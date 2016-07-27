/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 */
#include "ldlt_app.hxx"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <utility>

#include <omp.h>

#include "../AlignedAllocator.hxx"
#include "../BlockPool.hxx"
#include "../cpu_iface.hxx"
#include "block_ldlt.hxx"
#include "ldlt_tpp.hxx"
#include "common.hxx"
#include "wrappers.hxx"

namespace spral { namespace ssids { namespace cpu {

namespace ldlt_app_internal {

static const int INNER_BLOCK_SIZE = 32;

/** \return number of blocks for given n */
inline int calc_nblk(int n, int block_size) {
   return (n-1) / block_size + 1;
}

/** \return block size of block blk if maximum in dimension is n */
inline int calc_blkn(int blk, int n, int block_size) {
   return std::min(block_size, n-blk*block_size);
}

/** Workspace allocated on a per-thread basis */
template<typename T, typename Allocator=std::allocator<T>>
class ThreadWork {
   typedef typename std::allocator_traits<Allocator>::template rebind_traits<T> AllocTTraits;
   typedef typename std::allocator_traits<Allocator>::template rebind_traits<int> AllocIntTraits;
public:

   ThreadWork(ThreadWork const&) =delete;
   ThreadWork& operator=(ThreadWork const&) =delete;
   ThreadWork() =default;
   void init(int max_block_size, Allocator const& alloc = Allocator())
   {
      allocT_ = typename AllocTTraits::allocator_type(alloc);
      allocInt_ = typename AllocIntTraits::allocator_type(alloc);
      block_size_ = max_block_size;
      ld = AllocTTraits::allocate(allocT_, block_size_*block_size_);
      perm = AllocIntTraits::allocate(allocInt_, block_size_);
   }
   ~ThreadWork() {
      AllocIntTraits::deallocate(allocInt_, perm, block_size_);
      AllocTTraits::deallocate(allocT_, ld, block_size_*block_size_);
   }

public:
   T *ld;
   int *perm;

private:
   typename AllocTTraits::allocator_type allocT_;
   typename AllocIntTraits::allocator_type allocInt_;
   int block_size_;
};

template<typename T, typename IntAllocator>
class Column {
public:
   int nelim; //< Number of eliminated entries in this column
   T *d; //< pointer to local d
   std::vector<int, IntAllocator> lperm; //< local permutation

   Column(Column const&) {
      // Provided only for vector support, should never be used!
      omp_init_lock(&lock_);
   }
   Column& operator=(Column const&) =delete; // must be unique
   Column(int block_size, IntAllocator alloc=IntAllocator())
   : lperm(block_size, alloc)
   {
      omp_init_lock(&lock_);
   }
   ~Column() {
      omp_destroy_lock(&lock_);
   }

   /** Initialize number of passed columns ready for reduction */
   void init_passed(int passed) {
      npass_ = passed;
   }
   /** Updates number of passed columns (reduction by min) */
   void update_passed(int passed) {
      omp_set_lock(&lock_);
      npass_ = std::min(npass_, passed);
      omp_unset_lock(&lock_);
   }

   /** Adjust column after all blocks have passed to avoid split pivots */
   void adjust(int& next_elim) {
      // Test if last passed column was first part of a 2x2: if so,
      // decrement npass
      if(npass_>0) {
         T d11 = d[2*(npass_-1)+0];
         T d21 = d[2*(npass_-1)+1];
         if(std::isfinite(d11) && // not second half of 2x2
               d21 != 0.0)        // not a 1x1 or zero pivot
            npass_--;              // so must be first half 2x2
      }
      // Update elimination progress
      next_elim += npass_;
      nelim = npass_;
   }

   /** Moves perm for eliminated columns to elim_perm
    * (which may overlap from the front). Puts uneliminated variables in
    * failed_perm (no need for d with failed vars). */
   void move_back(int n, int* perm, int* elim_perm, int* failed_perm) {
      if(perm != elim_perm) { // Don't move if memory is identical
         for(int i=0; i<nelim; ++i)
            *(elim_perm++) = perm[i];
      }
      // Copy failed perm
      for(int i=nelim; i<n; ++i)
         *(failed_perm++) = perm[i];
   }

private:
   omp_lock_t lock_; //< Lock for altering npass
   int npass_; //< Reduction variable for nelim
};

/** Returns true if ptr is suitably aligned for AVX, false if not */
bool is_aligned(void* ptr) {
   const int align = 32;
   return (reinterpret_cast<uintptr_t>(ptr) % align == 0);
}

/** Move up eliminated entries to fill any gaps left by failed pivots
 *  within diagonal block.
 *  Note that out and aval may overlap. */
template<typename T, typename Column>
void move_up_diag(Column const& idata, Column const& jdata, T* out, T const* aval, int lda) {
   if(out == aval) return; // don't bother moving if memory is the same
   for(int j=0; j<jdata.nelim; ++j)
   for(int i=0; i<idata.nelim; ++i)
      out[j*lda+i] = aval[j*lda+i];
}

/** Move up eliminated entries to fill any gaps left by failed pivots
 *  within rectangular block of matrix.
 *  Note that out and aval may overlap. */
template<typename T, typename Column>
void move_up_rect(int m, int rfrom, Column const& jdata, T* out, T const* aval, int lda) {
   if(out == aval) return; // don't bother moving if memory is the same
   for(int j=0; j<jdata.nelim; ++j)
   for(int i=rfrom; i<m; ++i)
      out[j*lda+i] = aval[j*lda+i];
}

/** Copies failed rows and columns^T to specified locations */
template<typename T, typename Column>
void copy_failed_diag(int m, int n, Column const& idata, Column const& jdata, T* rout, T* cout, T* dout, int ldout, T const* aval, int lda) {
   /* copy rows */
   for(int j=0; j<jdata.nelim; ++j)
   for(int i=idata.nelim, iout=0; i<m; ++i, ++iout)
      rout[j*ldout+iout] = aval[j*lda+i];
   /* copy cols in transpose (not for diagonal block) */
   if(&idata != &jdata) {
      for(int j=jdata.nelim, jout=0; j<n; ++j, ++jout)
      for(int i=0; i<idata.nelim; ++i)
         cout[i*ldout+jout] = aval[j*lda+i];
   }
   /* copy intersection of failed rows and cols */
   for(int j=jdata.nelim, jout=0; j<n; j++, ++jout)
   for(int i=idata.nelim, iout=0; i<m; ++i, ++iout)
      dout[jout*ldout+iout] = aval[j*lda+i];
}

/** Copies failed columns to specified location */
template<typename T, typename Column>
void copy_failed_rect(int m, int n, int rfrom, Column const& jdata, T* cout, int ldout, T const* aval, int lda) {
   for(int j=jdata.nelim, jout=0; j<n; ++j, ++jout)
      for(int i=rfrom; i<m; ++i)
         cout[jout*ldout+i] = aval[j*lda+i];
}

/** Check if a block satisifies pivot threshold (colwise version) */
template <enum operation op, typename T>
int check_threshold(int rfrom, int rto, int cfrom, int cto, T u, T* aval, int lda) {
   // Perform threshold test for each uneliminated row/column
   for(int j=cfrom; j<cto; j++)
   for(int i=rfrom; i<rto; i++)
      if(fabs(aval[j*lda+i]) > 1.0/u)
         return (op==OP_N) ? j : i;
   // If we get this far, everything is good
   return (op==OP_N) ? cto : rto;
}

/** Performs solve with diagonal block \f$L_{21} = A_{21} L_{11}^{-T} D_1^{-1}\f$. Designed for below diagonal. */
/* NB: d stores (inverted) pivots as follows:
 * 2x2 ( a b ) stored as d = [ a b Inf c ]
 *     ( b c )
 * 1x1  ( a )  stored as d = [ a 0.0 ]
 * 1x1  ( 0 ) stored as d = [ 0.0 0.0 ]
 */
template <enum operation op, typename T>
void apply_pivot(int m, int n, int from, const T *diag, const T *d, const T small, T* aval, int lda) {
   if(op==OP_N && from > m) return; // no-op
   if(op==OP_T && from > n) return; // no-op

   if(op==OP_N) {
      // Perform solve L_11^-T
      host_trsm<T>(SIDE_RIGHT, FILL_MODE_LWR, OP_T, DIAG_UNIT,
            m, n, 1.0, diag, lda, aval, lda);
      // Perform solve L_21 D^-1
      for(int i=0; i<n; ) {
         if(i+1==n || std::isfinite(d[2*i+2])) {
            // 1x1 pivot
            T d11 = d[2*i];
            if(d11 == 0.0) {
               // Handle zero pivots carefully
               for(int j=0; j<m; j++) {
                  T v = aval[i*lda+j];
                  aval[i*lda+j] = 
                     (fabs(v)<small) ? 0.0
                                     : std::numeric_limits<T>::infinity()*v;
                  // NB: *v above handles NaNs correctly
               }
            } else {
               // Non-zero pivot, apply in normal fashion
               for(int j=0; j<m; j++)
                  aval[i*lda+j] *= d11;
            }
            i++;
         } else {
            // 2x2 pivot
            T d11 = d[2*i];
            T d21 = d[2*i+1];
            T d22 = d[2*i+3];
            for(int j=0; j<m; j++) {
               T a1 = aval[i*lda+j];
               T a2 = aval[(i+1)*lda+j];
               aval[i*lda+j]     = d11*a1 + d21*a2;
               aval[(i+1)*lda+j] = d21*a1 + d22*a2;
            }
            i += 2;
         }
      }
   } else { /* op==OP_T */
      // Perform solve L_11^-1
      host_trsm<T>(SIDE_LEFT, FILL_MODE_LWR, OP_N, DIAG_UNIT,
            m, n-from, 1.0, diag, lda, &aval[from*lda], lda);
      // Perform solve D^-T L_21^T
      for(int i=0; i<m; ) {
         if(i+1==m || std::isfinite(d[2*i+2])) {
            // 1x1 pivot
            T d11 = d[2*i];
            if(d11 == 0.0) {
               // Handle zero pivots carefully
               for(int j=from; j<n; j++) {
                  T v = aval[j*lda+i];
                  aval[j*lda+i] = 
                     (fabs(v)<small) ? 0.0 // *v handles NaNs
                                     : std::numeric_limits<T>::infinity()*v;
                  // NB: *v above handles NaNs correctly
               }
            } else {
               // Non-zero pivot, apply in normal fashion
               for(int j=from; j<n; j++) {
                  aval[j*lda+i] *= d11;
               }
            }
            i++;
         } else {
            // 2x2 pivot
            T d11 = d[2*i];
            T d21 = d[2*i+1];
            T d22 = d[2*i+3];
            for(int j=from; j<n; j++) {
               T a1 = aval[j*lda+i];
               T a2 = aval[j*lda+(i+1)];
               aval[j*lda+i]     = d11*a1 + d21*a2;
               aval[j*lda+(i+1)] = d21*a1 + d22*a2;
            }
            i += 2;
         }
      }
   }
}

/** Calculates LD from L and D */
template <enum operation op, typename T>
void calcLD(int m, int n, const T *l, int ldl, const T *d, T *ld, int ldld) {
   for(int col=0; col<n; ) {
      if(col+1==n || std::isfinite(d[2*col+2])) {
         // 1x1 pivot
         T d11 = d[2*col];
         if(d11 != 0.0) d11 = 1/d11; // Zero pivots just cause zeroes
         for(int row=0; row<m; row++)
            ld[col*ldld+row] = d11 * ((op==OP_N) ? l[col*ldl+row]
                                                 : l[row*ldl+col]);
         col++;
      } else {
         // 2x2 pivot
         T d11 = d[2*col];
         T d21 = d[2*col+1];
         T d22 = d[2*col+3];
         T det = d11*d22 - d21*d21;
         d11 = d11/det;
         d21 = d21/det;
         d22 = d22/det;
         for(int row=0; row<m; row++) {
            T a1 = (op==OP_N) ? l[col*ldl+row]     : l[row*ldl+col];
            T a2 = (op==OP_N) ? l[(col+1)*ldl+row] : l[row*ldl+(col+1)];
            ld[col*ldld+row]     =  d22*a1 - d21*a2;
            ld[(col+1)*ldld+row] = -d21*a1 + d11*a2;
         }
         col += 2;
      }
   }
}

template <typename T, typename Allocator=std::allocator<T*>>
class PoolBackup {
   typedef typename std::allocator_traits<Allocator>::template rebind_alloc<T*> TptrAlloc;
public:
   // FIXME: reduce pool size
   PoolBackup(int m, int n, int block_size, Allocator const& alloc=Allocator())
   : m_(m), n_(n), block_size_(block_size), mblk_(calc_nblk(m,block_size)),
     pool_(calc_nblk(n,block_size)*((calc_nblk(n,block_size)+1)/2+mblk_), block_size, alloc),
     ptr_(mblk_*calc_nblk(n,block_size), alloc)
   {}

   void acquire(int iblk, int jblk) {
      ptr_[jblk*mblk_+iblk] = pool_.get_wait();
   }

   void release(int iblk, int jblk) {
      pool_.release(ptr_[jblk*mblk_+iblk]);
      ptr_[jblk*mblk_+iblk] = nullptr;
   }

   void create_restore_point(int iblk, int jblk, T const* aval, int lda) {
      T*& lwork = ptr_[jblk*mblk_+iblk];
      lwork = pool_.get_wait();
      for(int j=0; j<get_ncol(jblk); j++)
      for(int i=0; i<get_nrow(iblk); i++)
         lwork[j*block_size_+i] = aval[j*lda+i];
   }

   /** Apply row permutation to block at same time as taking a copy */
   void create_restore_point_with_row_perm(int iblk, int jblk, int nperm, const int *lperm, T* aval, int lda) {
      T*& lwork = ptr_[jblk*mblk_+iblk];
      lwork = pool_.get_wait();
      for(int j=0; j<get_ncol(jblk); j++) {
         for(int i=0; i<nperm; i++) {
            int r = lperm[i];
            lwork[j*block_size_+i] = aval[j*lda+r];
         }
         for(int i=nperm; i<get_nrow(iblk); i++) {
            lwork[j*block_size_+i] = aval[j*lda+i];
         }
      }
      for(int j=0; j<get_ncol(jblk); j++)
      for(int i=0; i<nperm; i++)
         aval[j*lda+i] = lwork[j*block_size_+i];
   }

   /** Apply column permutation to block at same time as taking a copy */
   void create_restore_point_with_col_perm(int iblk, int jblk, const int *lperm, T* aval, int lda) {
      T*& lwork = ptr_[jblk*mblk_+iblk];
      lwork = pool_.get_wait();
      for(int j=0; j<get_ncol(jblk); j++) {
         int c = lperm[j];
         for(int i=0; i<get_nrow(iblk); i++)
            lwork[j*block_size_+i] = aval[c*lda+i];
      }
      for(int j=0; j<get_ncol(jblk); j++)
      for(int i=0; i<get_nrow(iblk); i++)
         aval[j*lda+i] = lwork[j*block_size_+i];
   }

   /** Restores any columns that have failed back to their previous
    *  values stored in lwork[] */
   void restore_part(int iblk, int jblk, int rfrom, int cfrom, T* aval, int lda) {
      T*& lwork = ptr_[jblk*mblk_+iblk];
      for(int j=cfrom; j<get_ncol(jblk); j++)
      for(int i=rfrom; i<get_nrow(iblk); i++)
         aval[j*lda+i] = lwork[j*block_size_+i];
   }

   /** Restores any columns that have failed back to their previous
    *  values stored in lwork[]. Applies a symmetric permutation while
    *  doing so. */
   void restore_part_with_sym_perm(int iblk, int jblk, int from, const int *lperm, T* aval, int lda) {
      T*& lwork = ptr_[jblk*mblk_+iblk];
      for(int j=from; j<get_ncol(jblk); j++) {
         int c = lperm[j];
         for(int i=from; i<get_ncol(jblk); i++) {
            int r = lperm[i];
            aval[j*lda+i] = (r>c) ? lwork[c*block_size_+r]
                                  : lwork[r*block_size_+c];
         }
         for(int i=get_ncol(jblk); i<get_nrow(iblk); i++)
            aval[j*lda+i] = lwork[c*block_size_+i];
      }
   }

private:
   inline int get_ncol(int blk) {
      return calc_blkn(blk, n_, block_size_);
   }
   inline int get_nrow(int blk) {
      return calc_blkn(blk, m_, block_size_);
   }

   int const m_;
   int const n_;
   int const block_size_;
   int const mblk_;
   BlockPool<T, Allocator> pool_;
   std::vector<T*, TptrAlloc> ptr_;
};

template<typename T,
         int BLOCK_SIZE,
         typename Backup,
         bool debug=false,
         typename Allocator=std::allocator<T>
         >
class LDLT;

template<typename T, int INNER_BLOCK_SIZE, typename ColAlloc>
class Block {
   typedef typename std::allocator_traits<ColAlloc>::template rebind_alloc<int> IntAlloc;
public:
   Block(int i, int j, int m, int n, std::vector<Column<T,IntAlloc>, ColAlloc>& cdata, T* a, int lda, int block_size)
   : i_(i), j_(j), m_(m), n_(n), lda_(lda), block_size_(block_size), cdata_(cdata),
     aval_(&a[j*block_size*lda+i*block_size])
   {}

   template <typename Backup>
   void backup(Backup& backup) {
      backup.create_restore_point(i_, j_, aval_, lda_);
   }

   template <typename Backup>
   void apply_rperm_and_backup(Backup& backup) {
      backup.create_restore_point_with_row_perm(
            i_, j_, get_ncol(i_), &cdata_[i_].lperm[0], aval_, lda_
            );
   }

   template <typename Backup>
   void apply_cperm_and_backup(Backup& backup) {
      backup.create_restore_point_with_col_perm(
            i_, j_, &cdata_[j_].lperm[0], aval_, lda_
            );
   }

   template <typename Backup>
   void restore_if_required(Backup& backup, int elim_col) {
      if(i_ == elim_col && j_ == elim_col) { // In eliminated diagonal block
         if(cdata_[i_].nelim < ncol()) { // If there are failed pivots
            backup.restore_part_with_sym_perm(
                  i_, j_, cdata_[i_].nelim, &cdata_[i_].lperm[0], aval_, lda_
                  );
         }
         // Release resources regardless, no longer required
         backup.release(i_, j_);
      }
      else if(i_ == elim_col) { // In eliminated row
         if(cdata_[i_].nelim < nrow()) // If there are failed pivots
            backup.restore_part(
                  i_, j_, cdata_[i_].nelim, cdata_[j_].nelim, aval_, lda_
                  );
         // Release resources regardless, no longer required
         backup.release(i_, j_);
      }
      else if(j_ == elim_col) { // In eliminated col
         if(cdata_[j_].nelim < ncol()) { // If there are failed pivots
            int rfrom = (i_ <= elim_col) ? cdata_[i_].nelim : 0;
            backup.restore_part(i_, j_, rfrom, cdata_[j_].nelim, aval_, lda_);
         }
         // Release resources regardless, no longer required
         backup.release(i_, j_);
      }
   }

   template <typename Allocator, typename ThreadWork>
   int factor(int& next_elim, int* perm, T* d, ThreadWork& work, struct cpu_factor_options const &options, Allocator const& alloc) {
      if(i_ != j_)
         throw std::runtime_error("factor called on non-diagonal block!");
      for(int i=0; i<ncol(); i++)
         cdata_[i_].lperm[i] = i;
      cdata_[i_].d = &d[2*next_elim];
      if(block_size_ != INNER_BLOCK_SIZE) {
         // Recurse
         PoolBackup<T, Allocator> inner_backup(
               nrow(), ncol(), INNER_BLOCK_SIZE, alloc
               );
         cdata_[i_].nelim =
            LDLT<T, INNER_BLOCK_SIZE, PoolBackup<T,Allocator>,
                 false, Allocator>
                :: factor(
                      nrow(), ncol(), &cdata_[i_].lperm[0], aval_, lda_,
                      cdata_[i_].d, inner_backup, options, INNER_BLOCK_SIZE,
                      alloc
                      );
         int *temp = work.perm;
         int* blkperm = &perm[i_*INNER_BLOCK_SIZE];
         for(int i=0; i<ncol(); ++i)
            temp[i] = blkperm[cdata_[i_].lperm[i]];
         for(int i=0; i<ncol(); ++i)
            blkperm[i] = temp[i];
         return cdata_[i_].nelim;
      } else { /* block_size == INNER_BLOCK_SIZE */
         // Call another routine for small block factorization
         if(ncol() < INNER_BLOCK_SIZE || !is_aligned(aval_)) {
            cdata_[i_].nelim = ldlt_tpp_factor(
                  nrow(), ncol(), &cdata_[i_].lperm[0], aval_, lda_,
                  cdata_[i_].d, work.ld, INNER_BLOCK_SIZE, options.u,
                  options.small
                  );
            int *temp = work.perm;
            int* blkperm = &perm[i_*INNER_BLOCK_SIZE];
            for(int i=0; i<ncol(); ++i)
               temp[i] = blkperm[cdata_[i_].lperm[i]];
            for(int i=0; i<ncol(); ++i)
               blkperm[i] = temp[i];
         } else {
            int* blkperm = &perm[i_*INNER_BLOCK_SIZE];
            block_ldlt<T, INNER_BLOCK_SIZE>(
                  0, blkperm, aval_, lda_, cdata_[i_].d, work.ld, options.u,
                  options.small, &cdata_[i_].lperm[0]
                  );
            cdata_[i_].nelim = INNER_BLOCK_SIZE;
         }
         return cdata_[i_].nelim;
      }
   }

   int apply_pivot_app(Block const& dblk, T u, T small) {
      if(i_ == j_)
         throw std::runtime_error("apply_pivot called on diagonal block!");
      if(i_ == dblk.i_) { // Apply within row (ApplyT)
         apply_pivot<OP_T>(
               cdata_[i_].nelim, ncol(), cdata_[j_].nelim, dblk.aval_,
               cdata_[i_].d, small, aval_, lda_
               );
         return check_threshold<OP_T>(
               0, cdata_[i_].nelim, cdata_[j_].nelim, ncol(), u, aval_, lda_
               );
      } else if(j_ == dblk.j_) { // Apply within column (ApplyN)
         apply_pivot<OP_N>(
               nrow(), cdata_[j_].nelim, 0, dblk.aval_,
               cdata_[j_].d, small, aval_, lda_
               );
         return check_threshold<OP_N>(
               0, nrow(), 0, cdata_[j_].nelim, u, aval_, lda_
               );
      } else {
         throw std::runtime_error("apply_pivot called on block outside eliminated column");
      }
   }

   template <typename ThreadWork>
   void update(Block const& isrc, Block const& jsrc, ThreadWork& work) {
      if(isrc.i_ == i_ && isrc.j_ == jsrc.j_) {
         // Update to right of elim column (UpdateN)
         int elim_col = isrc.j_;
         if(cdata_[elim_col].nelim == 0) return; // nothing to do
         int rfrom = (i_ <= elim_col) ? cdata_[i_].nelim : 0;
         int cfrom = (j_ <= elim_col) ? cdata_[j_].nelim : 0;
         calcLD<OP_N>(
               nrow()-rfrom, cdata_[elim_col].nelim, &isrc.aval_[rfrom],
               lda_, cdata_[elim_col].d, work.ld, block_size_
               );
         host_gemm(
               OP_N, OP_T, nrow()-rfrom, ncol()-cfrom,
               cdata_[elim_col].nelim, -1.0, work.ld, block_size_,
               &jsrc.aval_[cfrom], lda_, 1.0, &aval_[cfrom*lda_+rfrom], lda_
               );
      } else {
         // Update to left of elim column (UpdateT)
         int elim_col = jsrc.i_;
         if(cdata_[elim_col].nelim == 0) return; // nothing to do
         int rfrom = (i_ <= elim_col) ? cdata_[i_].nelim : 0;
         int cfrom = (j_ <= elim_col) ? cdata_[j_].nelim : 0;
         if(isrc.j_==elim_col) {
            calcLD<OP_N>(
                  nrow()-rfrom, cdata_[elim_col].nelim,
                  &isrc.aval_[rfrom], lda_,
                  cdata_[elim_col].d, work.ld, block_size_
                  );
         } else {
            calcLD<OP_T>(
                  nrow()-rfrom, cdata_[elim_col].nelim, &
                  isrc.aval_[rfrom*lda_], lda_,
                  cdata_[elim_col].d, work.ld, block_size_
                  );
         }
         host_gemm(
               OP_N, OP_N, nrow()-rfrom, ncol()-cfrom,
               cdata_[elim_col].nelim, -1.0, work.ld, block_size_,
               &jsrc.aval_[cfrom*lda_], lda_, 1.0, &aval_[cfrom*lda_+rfrom],
               lda_
               );
      }
   }

   int nrow() const { return get_nrow(i_); }
   int ncol() const { return get_ncol(j_); }
private:
   inline int get_ncol(int blk) const {
      return calc_blkn(blk, n_, block_size_);
   }
   inline int get_nrow(int blk) const {
      return calc_blkn(blk, m_, block_size_);
   }

   int const i_; //< block's row
   int const j_; //< block's column
   int const m_; //< global number of rows
   int const n_; //< global number of columns
   int const lda_; //< leading dimension of underlying storage
   int const block_size_; //< block size
   std::vector<Column<T,IntAlloc>, ColAlloc>& cdata_; //< global column data array
   T* aval_;
};

template<typename T,
         int BLOCK_SIZE,
         typename Backup,
         bool debug,
         typename Allocator
         >
class LDLT {
   typedef typename std::allocator_traits<Allocator>::template rebind_alloc<int> IntAlloc;
   typedef typename std::allocator_traits<Allocator>::template rebind_alloc<Column<T,IntAlloc>> ColAlloc;
   typedef typename std::allocator_traits<Allocator>::template rebind_alloc<bool> BoolAlloc;
   typedef typename std::allocator_traits<Allocator>::template rebind_alloc<T> TAlloc;
private:
   template <typename ThreadWork>
   static
   int run_elim(int const m, int const n, int* perm, T* a, int const lda, T* d,
         std::vector<Column<T,IntAlloc>, ColAlloc>& cdata, Backup& backup,
         ThreadWork all_thread_work[],
         struct cpu_factor_options const& options, int const block_size,
         Allocator const& alloc) {
      typedef Block<T, BLOCK_SIZE, ColAlloc> BlockSpec;
      //printf("ENTRY %d %d vis %d %d %d\n", m, n, mblk, nblk, block_size);

      int const nblk = calc_nblk(n, block_size);
      int const mblk = calc_nblk(m, block_size);

      /* Setup */
      int next_elim = 0;

      /* Inner loop - iterate over block columns */
      for(int blk=0; blk<nblk; blk++) {
         /*if(debug) {
            printf("Bcol %d:\n", blk);
            print_mat(mblk, nblk, m, n, blkdata, cdata, lda);
         }*/

         // Factor diagonal: depend on perm[blk*block_size] as we init npass
         #pragma omp task default(none) \
            firstprivate(blk) \
            shared(a, perm, backup, cdata, all_thread_work, next_elim, d, \
                   options, alloc) \
            depend(inout: a[blk*block_size*lda+blk*block_size:1]) \
            depend(inout: perm[blk*block_size:1])
         {
            if(debug) printf("Factor(%d)\n", blk);
            int thread_num = omp_get_thread_num();
            BlockSpec dblk(blk, blk, m, n, cdata, a, lda, block_size);
            // Store a copy for recovery in case of a failed column
            dblk.backup(backup);
            // Perform actual factorization
            int nelim = dblk.template factor<Allocator>(
                  next_elim, perm, d, all_thread_work[thread_num],
                  options, alloc
                  );
            // Init threshold check (non locking => task dependencies)
            cdata[blk].init_passed(nelim);
         }
         
         // Loop over off-diagonal blocks applying pivot
         for(int jblk=0; jblk<blk; jblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, jblk) \
               shared(a, backup, cdata, options) \
               depend(in: a[blk*block_size*lda+blk*block_size:1]) \
               depend(inout: a[jblk*block_size*lda+blk*block_size:1]) \
               depend(in: perm[blk*block_size:1])
            {
               if(debug) printf("ApplyT(%d,%d)\n", blk, jblk);
               BlockSpec dblk(blk, blk, m, n, cdata, a, lda, block_size);
               BlockSpec cblk(blk, jblk, m, n, cdata, a, lda, block_size);
               // Apply row permutation from factorization of dblk and in
               // the process, store a (permuted) copy for recovery in case of
               // a failed column
               cblk.apply_rperm_and_backup(backup);
               // Perform elimination and determine number of rows in block
               // passing a posteori threshold pivot test
               int blkpass = cblk.apply_pivot_app(
                     dblk, options.u, options.small
                     );
               // Update column's passed pivot count
               cdata[blk].update_passed(blkpass);
            }
         }
         for(int iblk=blk+1; iblk<mblk; iblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, iblk) \
               shared(a, backup, cdata, options) \
               depend(in: a[blk*block_size*lda+blk*block_size:1]) \
               depend(inout: a[blk*block_size*lda+iblk*block_size:1]) \
               depend(in: perm[blk*block_size:1])
            {
               if(debug) printf("ApplyN(%d,%d)\n", iblk, blk);
               BlockSpec dblk(blk, blk, m, n, cdata, a, lda, block_size);
               BlockSpec rblk(iblk, blk, m, n, cdata, a, lda, block_size);
               // Apply column permutation from factorization of dblk and in
               // the process, store a (permuted) copy for recovery in case of
               // a failed column
               rblk.apply_cperm_and_backup(backup);
               // Perform elimination and determine number of rows in block
               // passing a posteori threshold pivot test
               int blkpass = rblk.apply_pivot_app(dblk, options.u, options.small);
               // Update column's passed pivot count
               cdata[blk].update_passed(blkpass);
            }
         }

         // Adjust column once all applys have finished and we know final
         // number of passed columns.
         #pragma omp task default(none) \
            firstprivate(blk) \
            shared(cdata, next_elim) \
            depend(inout: perm[blk*block_size:1])
         {
            if(debug) printf("Adjust(%d)\n", blk);
            cdata[blk].adjust(next_elim);
         }

         // Update uneliminated columns
         for(int jblk=0; jblk<blk; jblk++) {
            for(int iblk=jblk; iblk<mblk; iblk++) {
               // Calculate block index we depend on for i
               // (we only work with lower half of matrix)
               int adep_idx = (blk<iblk) ? blk*block_size*lda + iblk*block_size
                                         : iblk*block_size*lda + blk*block_size;
               #pragma omp task default(none) \
                  firstprivate(blk, iblk, jblk) \
                  shared(a, cdata, backup, all_thread_work)\
                  depend(inout: a[jblk*block_size*lda+iblk*block_size:1]) \
                  depend(in: perm[blk*block_size:1]) \
                  depend(in: a[jblk*block_size*lda+blk*block_size:1]) \
                  depend(in: a[adep_idx:1])
               {
                  if(debug) printf("UpdateT(%d,%d,%d)\n", iblk, jblk, blk);
                  int thread_num = omp_get_thread_num();
                  BlockSpec ublk(iblk, jblk, m, n, cdata, a, lda, block_size);
                  int isrc_row = (blk<=iblk) ? iblk : blk;
                  int isrc_col = (blk<=iblk) ? blk : iblk;
                  BlockSpec isrc(isrc_row, isrc_col, m, n, cdata, a, lda,
                        block_size);
                  BlockSpec jsrc(blk, jblk, m, n, cdata, a, lda, block_size);
                  // If we're on the block row we've just eliminated, restore
                  // any failed rows and release resources storing backup
                  ublk.restore_if_required(backup, blk);
                  // Perform actual update
                  ublk.update(isrc, jsrc, all_thread_work[thread_num]);
               }
            }
         }
         for(int jblk=blk; jblk<nblk; jblk++) {
            for(int iblk=jblk; iblk<mblk; iblk++) {
               #pragma omp task default(none) \
                  firstprivate(blk, iblk, jblk) \
                  shared(a, cdata, backup, all_thread_work) \
                  depend(inout: a[jblk*block_size*lda+iblk*block_size:1]) \
                  depend(in: perm[blk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+iblk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+jblk*block_size:1])
               {
                  if(debug) printf("UpdateN(%d,%d,%d)\n", iblk, jblk, blk);
                  int thread_num = omp_get_thread_num();
                  BlockSpec ublk(iblk, jblk, m, n, cdata, a, lda, block_size);
                  BlockSpec isrc(iblk, blk, m, n, cdata, a, lda, block_size);
                  BlockSpec jsrc(jblk, blk, m, n, cdata, a, lda, block_size);
                  // If we're on the block col we've just eliminated, restore
                  // any failed cols and release resources storing backup
                  ublk.restore_if_required(backup, blk);
                  // Perform actual update
                  ublk.update(isrc, jsrc, all_thread_work[thread_num]);
               }
            }
         }
      }
      #pragma omp taskwait

      /*if(debug) {
         printf("PostElim:\n");
         print_mat(mblk, nblk, m, n, blkdata, cdata, lda);
      }*/

      return next_elim;
   }

   static
   void print_mat(int m, int n, const int *perm, std::vector<bool, BoolAlloc> const& eliminated, const T *a, int lda) {
      for(int row=0; row<m; row++) {
         if(row < n)
            printf("%d%s:", perm[row], eliminated[row]?"X":" ");
         else
            printf("%d%s:", row, "U");
         for(int col=0; col<std::min(n,row+1); col++)
            printf(" %10.4f", a[col*lda+row]);
         printf("\n");
      }
   }

   static
   inline int get_ncol(int blk, int n, int block_size) {
      return calc_blkn(blk, n, block_size);
   }
   static
   inline int get_nrow(int blk, int m, int block_size) {
      return calc_blkn(blk, m, block_size);
   }

public:
   /** Factorize an entire matrix */
   static
   int factor(int m, int n, int *perm, T *a, int lda, T *d, Backup& backup, struct cpu_factor_options const& options, int block_size, Allocator const& alloc=Allocator()) {
      /* Sanity check arguments */
      if(m < n) return -1;
      if(lda < n) return -4;

      /* Initialize useful quantities: */
      int nblk = calc_nblk(n, block_size);
      int mblk = calc_nblk(m, block_size);

      /* Temporary workspaces */
      int num_threads = omp_get_max_threads();
      ThreadWork<T,Allocator> all_thread_work[num_threads];
      for(int i=0; i<num_threads; ++i)
         all_thread_work[i].init(block_size, alloc);
      std::vector<Column<T,IntAlloc>, ColAlloc> cdata(alloc);
      cdata.reserve(nblk);
      for(int i=0; i<nblk; ++i)
         cdata.emplace_back(block_size, alloc);

      /* Main loop
       *    - Each pass leaves any failed pivots in place and keeps everything
       *      up-to-date.
       *    - If no pivots selected across matrix, perform swaps to get large
       *      entries into diagonal blocks
       */
      int num_elim = run_elim(
            m, n, perm, a, lda, d, cdata, backup, all_thread_work, options,
            block_size, alloc
            );

      // Permute failed entries to end
      std::vector<int, IntAlloc> failed_perm(n-num_elim, alloc);
      for(int jblk=0, insert=0, fail_insert=0; jblk<nblk; jblk++) {
         cdata[jblk].move_back(
               get_ncol(jblk, n, block_size), &perm[jblk*block_size],
               &perm[insert], &failed_perm[fail_insert]
               );
         insert += cdata[jblk].nelim;
         fail_insert += get_ncol(jblk, n, block_size) - cdata[jblk].nelim;
      }
      for(int i=0; i<n-num_elim; ++i)
         perm[num_elim+i] = failed_perm[i];

      // Extract failed entries of a
      int nfail = n-num_elim;
      std::vector<T, TAlloc> failed_diag(nfail*n, alloc);
      std::vector<T, TAlloc> failed_rect(nfail*(m-n), alloc);
      for(int jblk=0, jfail=0, jinsert=0; jblk<nblk; ++jblk) {
         // Diagonal part
         for(int iblk=jblk, ifail=jfail, iinsert=jinsert; iblk<nblk; ++iblk) {
            copy_failed_diag(
                  get_ncol(iblk, n, block_size), get_ncol(jblk, n, block_size),
                  cdata[iblk], cdata[jblk],
                  &failed_diag[jinsert*nfail+ifail],
                  &failed_diag[iinsert*nfail+jfail],
                  &failed_diag[num_elim*nfail+jfail*nfail+ifail],
                  nfail, &a[jblk*block_size*lda+iblk*block_size], lda
                  );
            iinsert += cdata[iblk].nelim;
            ifail += get_ncol(iblk, n, block_size) - cdata[iblk].nelim;
         }
         // Rectangular part
         // (be careful with blocks that contain both diag and rect parts)
         copy_failed_rect(
               get_nrow(nblk-1, m, block_size), get_ncol(jblk, n, block_size),
               get_ncol(nblk-1, n, block_size), cdata[jblk],
               &failed_rect[jfail*(m-n)+(nblk-1)*block_size-n], m-n,
               &a[jblk*block_size*lda+(nblk-1)*block_size], lda
               );
         for(int iblk=nblk; iblk<mblk; ++iblk) {
            copy_failed_rect(
                  get_nrow(iblk, m, block_size),
                  get_ncol(jblk, n, block_size), 0, cdata[jblk],
                  &failed_rect[jfail*(m-n)+iblk*block_size-n], m-n,
                  &a[jblk*block_size*lda+iblk*block_size], lda
                  );
         }
         jinsert += cdata[jblk].nelim;
         jfail += get_ncol(jblk, n, block_size) - cdata[jblk].nelim;
      }

      // Move data up
      for(int jblk=0, jinsert=0; jblk<nblk; ++jblk) {
         // Diagonal part
         for(int iblk=jblk, iinsert=jinsert; iblk<nblk; ++iblk) {
            move_up_diag(
                  cdata[iblk], cdata[jblk], &a[jinsert*lda+iinsert],
                  &a[jblk*block_size*lda+iblk*block_size], lda
                  );
            iinsert += cdata[iblk].nelim;
         }
         // Rectangular part
         // (be careful with blocks that contain both diag and rect parts)
         move_up_rect(
               get_nrow(nblk-1, m, block_size),
               get_ncol(nblk-1, n, block_size), cdata[jblk],
               &a[jinsert*lda+(nblk-1)*block_size],
               &a[jblk*block_size*lda+(nblk-1)*block_size], lda
               );
         for(int iblk=nblk; iblk<mblk; ++iblk)
            move_up_rect(
                  get_nrow(iblk, m, block_size), 0, cdata[jblk],
                  &a[jinsert*lda+iblk*block_size],
                  &a[jblk*block_size*lda+iblk*block_size], lda
                  );
         jinsert += cdata[jblk].nelim;
      }
      
      // Store failed entries back to correct locations
      // Diagonal part
      for(int j=0; j<n; ++j)
      for(int i=std::max(j,num_elim), k=i-num_elim; i<n; ++i, ++k)
         a[j*lda+i] = failed_diag[j*nfail+k];
      // Rectangular part
      T* arect = &a[num_elim*lda+n];
      for(int j=0; j<nfail; ++j)
      for(int i=0; i<m-n; ++i)
         arect[j*lda+i] = failed_rect[j*(m-n)+i];

      if(debug) {
         std::vector<bool, BoolAlloc> eliminated(n, alloc);
         for(int i=0; i<num_elim; i++) eliminated[i] = true;
         for(int i=num_elim; i<n; i++) eliminated[i] = false;
         printf("FINAL:\n");
         print_mat(m, n, perm, eliminated, a, lda);
      }

      return num_elim;
   }
};

} /* namespace spral::ssids:cpu::ldlt_app_internal */

using namespace spral::ssids::cpu::ldlt_app_internal;

template<typename T>
int ldlt_app_factor(int m, int n, int *perm, T *a, int lda, T *d, struct cpu_factor_options const& options) {
   // Things will break if outer block size is not a multiple of inner one
   if(options.cpu_task_block_size % INNER_BLOCK_SIZE != 0)
      throw std::runtime_error("options.cpu_task_block_size must be multiple of inner block size");

   // Template parameters and workspaces
   int const outer_block_size = INNER_BLOCK_SIZE;//options.cpu_task_block_size;
   bool const debug = false;
   typedef std::allocator<T> Allocator;
   Allocator alloc;
   PoolBackup<T, Allocator> backup(m, n, outer_block_size, alloc);

   // Actual call
   return LDLT
      <T, INNER_BLOCK_SIZE, PoolBackup<T,Allocator>, debug,
       Allocator>
      ::factor(
            m, n, perm, a, lda, d, backup, options,
            outer_block_size, alloc
            );
}
template int ldlt_app_factor<double>(int, int, int*, double*, int, double*, struct cpu_factor_options const&);

template <typename T>
void ldlt_app_solve_fwd(int m, int n, T const* l, int ldl, int nrhs, T* x, int ldx) {
   if(nrhs==1) {
      host_trsv(FILL_MODE_LWR, OP_N, DIAG_UNIT, n, l, ldl, x, 1);
      if(m > n)
         gemv(OP_N, m-n, n, -1.0, &l[n], ldl, x, 1, 1.0, &x[n], 1);
   } else {
      host_trsm(SIDE_LEFT, FILL_MODE_LWR, OP_N, DIAG_UNIT, n, nrhs, 1.0, l, ldl, x, ldx);
      if(m > n)
         host_gemm(OP_N, OP_N, m-n, nrhs, n, -1.0, &l[n], ldl, x, ldx, 1.0, &x[n], ldx);
   }
}
template void ldlt_app_solve_fwd<double>(int, int, double const*, int, int, double*, int);

template <typename T>
void ldlt_app_solve_diag(int n, T const* d, T* x) {
   for(int i=0; i<n; ) {
      if(i+1==n || std::isfinite(d[2*i+2])) {
         // 1x1 pivot
         T d11 = d[2*i];
         x[i] *= d11;
         i++;
      } else {
         // 2x2 pivot
         T d11 = d[2*i];
         T d21 = d[2*i+1];
         T d22 = d[2*i+3];
         T x1 = x[i];
         T x2 = x[i+1];
         x[i]   = d11*x1 + d21*x2;
         x[i+1] = d21*x1 + d22*x2;
         i += 2;
      }
   }
}
template void ldlt_app_solve_diag<double>(int, double const*, double*);

template <typename T>
void ldlt_app_solve_bwd(int m, int n, T const* l, int ldl, int nrhs, T* x, int ldx) {
   if(nrhs==1) {
      if(m > n)
         gemv(OP_T, m-n, n, -1.0, &l[n], ldl, &x[n], 1, 1.0, x, 1);
      host_trsv(FILL_MODE_LWR, OP_T, DIAG_UNIT, n, l, ldl, x, 1);
   } else {
      if(m > n)
         host_gemm(OP_T, OP_N, n, nrhs, m-n, -1.0, &l[n], ldl, &x[n], ldx, 1.0, x, ldx);
      host_trsm(SIDE_LEFT, FILL_MODE_LWR, OP_T, DIAG_UNIT, n, nrhs, 1.0, l, ldl, x, ldx);
   }
}
template void ldlt_app_solve_bwd<double>(int, int, double const*, int, int, double*, int);

}}} /* namespaces spral::ssids::cpu */
