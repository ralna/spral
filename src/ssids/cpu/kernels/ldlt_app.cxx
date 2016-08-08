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
#include <memory>
#include <ostream>
#include <sstream>
#include <utility>

#include <omp.h>

#include "../AlignedAllocator.hxx"
#include "../BlockPool.hxx"
#include "../BuddyAllocator.hxx"
#include "../cpu_iface.hxx"
#include "../Workspace.hxx"
#include "block_ldlt.hxx"
#include "calc_ld.hxx"
#include "ldlt_tpp.hxx"
#include "common.hxx"
#include "wrappers.hxx"

#include "../profile.hxx"

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

template<typename T>
class Column {
public:
   bool first_elim; //< True if first column with eliminations
   int nelim; //< Number of eliminated entries in this column
   T *d; //< pointer to local d

   Column(Column const&) =delete; // must be unique
   Column& operator=(Column const&) =delete; // must be unique
   Column() =default;

   /** Initialize number of passed columns ready for reduction */
   void init_passed(int passed) {
      npass_ = passed;
   }
   /** Updates number of passed columns (reduction by min) */
   void update_passed(int passed) {
      lock_.set();
      npass_ = std::min(npass_, passed);
      lock_.unset();
   }
   /** Return true if passed < nelim */
   bool test_fail(int passed) {
      bool fail = (passed < nelim);
      if(!fail) {
         // Record number of blocks in column passing this test
         #pragma omp atomic update
         ++npass_;
      }
      return fail;
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
      first_elim = (next_elim==0 && npass_>0);
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

   int get_npass() const { return npass_; }

private:
   spral::omp::Lock lock_; //< Lock for altering npass
   int npass_=0; //< Reduction variable for nelim
};

template<typename T, typename IntAlloc>
class ColumnData {
   typedef typename std::allocator_traits<IntAlloc>::template rebind_traits<Column<T>> ColAllocTraits;
   typedef typename std::allocator_traits<IntAlloc> IntAllocTraits;
public:
   ColumnData(ColumnData const&) =delete; //not copyable
   ColumnData& operator=(ColumnData const&) =delete; //not copyable
   ColumnData(int n, int block_size, IntAlloc const& alloc)
   : n_(n), block_size_(block_size), alloc_(alloc)
   {
      int nblk = calc_nblk(n_, block_size_);
      typename ColAllocTraits::allocator_type colAlloc(alloc_);
      cdata_ = ColAllocTraits::allocate(colAlloc, nblk);
      for(int i=0; i<nblk; ++i)
         ColAllocTraits::construct(colAlloc, &cdata_[i]);
      lperm_ = IntAllocTraits::allocate(alloc_, nblk*block_size_);
   }
   ~ColumnData() {
      int nblk = calc_nblk(n_, block_size_);
      IntAllocTraits::deallocate(alloc_, lperm_, nblk*block_size_);
      typename ColAllocTraits::allocator_type colAlloc(alloc_);
      ColAllocTraits::deallocate(colAlloc, cdata_, nblk);
   }

   Column<T>& operator[](int idx) { return cdata_[idx]; }

   int* get_lperm(int blk) { return &lperm_[blk*block_size_]; }

   /** Calculate number of eliminated columns in unpivoted case */
   int calc_nelim(int m) const {
      int mblk = calc_nblk(m, block_size_);
      int nblk = calc_nblk(n_, block_size_);
      int nelim = 0;
      for(int j=0; j<nblk; ++j) {
         if(cdata_[j].get_npass() == mblk-j)
            nelim += cdata_[j].nelim;
      }
      return nelim;
   };

private:
   int const n_;
   int const block_size_;
   IntAlloc alloc_;
   Column<T> *cdata_;
   int* lperm_;
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
   int least_fail = (op==OP_N) ? cto : rto;
   for(int j=cfrom; j<cto; j++)
   for(int i=rfrom; i<rto; i++)
      if(fabs(aval[j*lda+i]) > 1.0/u) {
         if(op==OP_N) {
            // must be least failed col
            return j;
         } else {
            // may be an earlier failed row
            least_fail = std::min(least_fail, i);
            break;
         }
      }
   // If we get this far, everything is good
   return least_fail;
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

template <typename T, typename Allocator=std::allocator<T>>
class CopyBackup {
public:
   CopyBackup(CopyBackup const&) =delete;
   CopyBackup& operator=(CopyBackup const&) =delete;
   CopyBackup(int m, int n, int block_size, Allocator const& alloc=Allocator())
   : alloc_(alloc), m_(m), n_(n), block_size_(block_size),
     ldcopy_(align_lda<T>(m_)), acopy_(alloc_.allocate(n_*ldcopy_))
   {}
   ~CopyBackup() {
      release_all_memory();
   }

   /** Release all underlying memory - instnace cannot be used again */
   void release_all_memory() {
      if(acopy_) {
         alloc_.deallocate(acopy_, n_*ldcopy_);
         acopy_ = nullptr;
      }
   }

   void release(int iblk, int jblk) { /* no-op */ }

   void create_restore_point(int iblk, int jblk, T const* aval, int lda) {
      T* lwork = get_lwork(iblk, jblk);
      for(int j=0; j<get_ncol(jblk); j++)
      for(int i=0; i<get_nrow(iblk); i++)
         lwork[j*ldcopy_+i] = aval[j*lda+i];
   }

   /** Apply row permutation to block at same time as taking a copy */
   void create_restore_point_with_row_perm(int iblk, int jblk, int nperm, const int *lperm, T* aval, int lda) {
      T* lwork = get_lwork(iblk, jblk);
      for(int j=0; j<get_ncol(jblk); j++) {
         for(int i=0; i<nperm; i++) {
            int r = lperm[i];
            lwork[j*ldcopy_+i] = aval[j*lda+r];
         }
         for(int i=nperm; i<get_nrow(iblk); i++) {
            lwork[j*ldcopy_+i] = aval[j*lda+i];
         }
      }
      for(int j=0; j<get_ncol(jblk); j++)
      for(int i=0; i<nperm; i++)
         aval[j*lda+i] = lwork[j*ldcopy_+i];
   }

   /** Apply column permutation to block at same time as taking a copy */
   void create_restore_point_with_col_perm(int iblk, int jblk, const int *lperm, T* aval, int lda) {
      T* lwork = get_lwork(iblk, jblk);
      for(int j=0; j<get_ncol(jblk); j++) {
         int c = lperm[j];
         for(int i=0; i<get_nrow(iblk); i++)
            lwork[j*ldcopy_+i] = aval[c*lda+i];
      }
      for(int j=0; j<get_ncol(jblk); j++)
      for(int i=0; i<get_nrow(iblk); i++)
         aval[j*lda+i] = lwork[j*ldcopy_+i];
   }

   /** Restores any columns that have failed back to their previous
    *  values stored in lwork[] */
   void restore_part(int iblk, int jblk, int rfrom, int cfrom, T* aval, int lda) {
      T* lwork = get_lwork(iblk, jblk);
      for(int j=cfrom; j<get_ncol(jblk); j++)
      for(int i=rfrom; i<get_nrow(iblk); i++)
         aval[j*lda+i] = lwork[j*ldcopy_+i];
   }

   /** Restores any columns that have failed back to their previous
    *  values stored in lwork[]. Applies a symmetric permutation while
    *  doing so. */
   void restore_part_with_sym_perm(int iblk, int jblk, int from, const int *lperm, T* aval, int lda) {
      T* lwork = get_lwork(iblk, jblk);
      for(int j=from; j<get_ncol(jblk); j++) {
         int c = lperm[j];
         for(int i=from; i<get_ncol(jblk); i++) {
            int r = lperm[i];
            aval[j*lda+i] = (r>c) ? lwork[c*ldcopy_+r]
                                  : lwork[r*ldcopy_+c];
         }
         for(int i=get_ncol(jblk); i<get_nrow(iblk); i++)
            aval[j*lda+i] = lwork[c*ldcopy_+i];
      }
   }

private:
   inline T* get_lwork(int iblk, int jblk) {
      return &acopy_[jblk*block_size_*ldcopy_+iblk*block_size_];
   }
   inline int get_ncol(int blk) const {
      return calc_blkn(blk, n_, block_size_);
   }
   inline int get_nrow(int blk) const {
      return calc_blkn(blk, m_, block_size_);
   }

   Allocator alloc_;
   int const m_;
   int const n_;
   int const block_size_;
   size_t const ldcopy_;
   T *acopy_;
};

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
         bool use_tasks, // Use tasks, so we can disable on one or more levels
         bool debug=false,
         typename Allocator=std::allocator<T>
         >
class LDLT;

template<typename T, int INNER_BLOCK_SIZE, typename IntAlloc>
class Block {
public:
   Block(int i, int j, int m, int n, ColumnData<T,IntAlloc>& cdata, T* a, int lda, int block_size)
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
            i_, j_, get_ncol(i_), cdata_.get_lperm(i_), aval_, lda_
            );
   }

   void apply_rperm(Workspace& work) {
      int ldl = align_lda<T>(block_size_);
      T* lwork = work.get_ptr<T>(ncol()*ldl);
      int* lperm = cdata_.get_lperm(i_);
      // Copy into lwork with permutation
      for(int j=0; j<ncol(); ++j) {
         for(int i=0; i<get_ncol(i_); ++i) {
            int r = lperm[i];
            lwork[j*ldl+i] = aval_[j*lda_+r];
         }
      }
      // Copy back again
      for(int j=0; j<ncol(); ++j)
      for(int i=0; i<get_ncol(i_); ++i)
         aval_[j*lda_+i] = lwork[j*ldl+i];
   }

   template <typename Backup>
   void apply_cperm_and_backup(Backup& backup) {
      backup.create_restore_point_with_col_perm(
            i_, j_, cdata_.get_lperm(j_), aval_, lda_
            );
   }

   void apply_cperm(Workspace& work) {
      int ldl = align_lda<T>(block_size_);
      T* lwork = work.get_ptr<T>(ncol()*ldl);
      int* lperm = cdata_.get_lperm(j_);
      // Copy into lwork with permutation
      for(int j=0; j<ncol(); ++j) {
         int c = lperm[j];
         for(int i=0; i<nrow(); ++i)
            lwork[j*ldl+i] = aval_[c*lda_+i];
      }
      // Copy back again
      for(int j=0; j<ncol(); ++j)
      for(int i=0; i<nrow(); ++i)
         aval_[j*lda_+i] = lwork[j*ldl+i];
   }

   template <typename Backup>
   void full_restore(Backup& backup) {
      backup.restore_part(i_, j_, 0, 0, aval_, lda_);
   }

   template <typename Backup>
   void restore_if_required(Backup& backup, int elim_col) {
      if(i_ == elim_col && j_ == elim_col) { // In eliminated diagonal block
         if(cdata_[i_].nelim < ncol()) { // If there are failed pivots
            backup.restore_part_with_sym_perm(
                  i_, j_, cdata_[i_].nelim, cdata_.get_lperm(i_), aval_, lda_
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

   template <typename Allocator>
   int factor(int& next_elim, int* perm, T* d, struct cpu_factor_options const &options, std::vector<Workspace>& work, Allocator const& alloc) {
      if(i_ != j_)
         throw std::runtime_error("factor called on non-diagonal block!");
      int* lperm = cdata_.get_lperm(i_);
      for(int i=0; i<ncol(); i++)
         lperm[i] = i;
      cdata_[i_].d = &d[2*next_elim];
      if(block_size_ != INNER_BLOCK_SIZE) {
         // Recurse
         CopyBackup<T, Allocator> inner_backup(
               nrow(), ncol(), INNER_BLOCK_SIZE, alloc
               );
         bool const use_tasks = false; // Don't run in parallel at lower level
         bool const debug = false; // Don't print debug info for inner call
         cdata_[i_].nelim =
            LDLT<T, INNER_BLOCK_SIZE, CopyBackup<T,Allocator>,
                 use_tasks, debug, Allocator>
                ::factor(
                      nrow(), ncol(), lperm, aval_, lda_,
                      cdata_[i_].d, inner_backup, options, options.pivot_method,
                      INNER_BLOCK_SIZE, 0, nullptr, 0, work, alloc
                      );
         int* temp = work[omp_get_thread_num()].get_ptr<int>(ncol());
         int* blkperm = &perm[i_*block_size_];
         for(int i=0; i<ncol(); ++i)
            temp[i] = blkperm[lperm[i]];
         for(int i=0; i<ncol(); ++i)
            blkperm[i] = temp[i];
      } else { /* block_size == INNER_BLOCK_SIZE */
         // Call another routine for small block factorization
         if(ncol() < INNER_BLOCK_SIZE || !is_aligned(aval_)) {
            T* ld = work[omp_get_thread_num()].get_ptr<T>(2*INNER_BLOCK_SIZE);
            cdata_[i_].nelim = ldlt_tpp_factor(
                  nrow(), ncol(), lperm, aval_, lda_,
                  cdata_[i_].d, ld, INNER_BLOCK_SIZE, options.u,
                  options.small
                  );
            int* temp = work[omp_get_thread_num()].get_ptr<int>(ncol());
            int* blkperm = &perm[i_*INNER_BLOCK_SIZE];
            for(int i=0; i<ncol(); ++i)
               temp[i] = blkperm[lperm[i]];
            for(int i=0; i<ncol(); ++i)
               blkperm[i] = temp[i];
         } else {
            int* blkperm = &perm[i_*INNER_BLOCK_SIZE];
            T* ld = work[omp_get_thread_num()].get_ptr<T>(
                  INNER_BLOCK_SIZE*INNER_BLOCK_SIZE
                  );
            block_ldlt<T, INNER_BLOCK_SIZE>(
                  0, blkperm, aval_, lda_, cdata_[i_].d, ld, options.u,
                  options.small, lperm
                  );
            cdata_[i_].nelim = INNER_BLOCK_SIZE;
         }
      }
      return cdata_[i_].nelim;
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

   void update(Block const& isrc, Block const& jsrc, Workspace& work, double beta=1.0, T* upd=nullptr, int ldupd=0) {
      if(isrc.i_ == i_ && isrc.j_ == jsrc.j_) {
         // Update to right of elim column (UpdateN)
         int elim_col = isrc.j_;
         if(cdata_[elim_col].nelim == 0) return; // nothing to do
         int rfrom = (i_ <= elim_col) ? cdata_[i_].nelim : 0;
         int cfrom = (j_ <= elim_col) ? cdata_[j_].nelim : 0;
         int ldld = align_lda<T>(block_size_);
         T* ld = work.get_ptr<T>(block_size_*ldld);
         // NB: we use ld[rfrom] below so alignment matches that of aval[rfrom]
         calcLD<OP_N>(
               nrow()-rfrom, cdata_[elim_col].nelim, &isrc.aval_[rfrom],
               lda_, cdata_[elim_col].d, &ld[rfrom], ldld
               );
         host_gemm(
               OP_N, OP_T, nrow()-rfrom, ncol()-cfrom, cdata_[elim_col].nelim,
               -1.0, &ld[rfrom], ldld, &jsrc.aval_[cfrom], lda_,
               1.0, &aval_[cfrom*lda_+rfrom], lda_
               );
         if(upd && j_==calc_nblk(n_,block_size_)-1) {
            // Handle fractional part of upd that "belongs" to this block
            int u_ncol = std::min(block_size_-ncol(), m_-n_); // ncol for upd
            beta = (cdata_[elim_col].first_elim) ? beta : 1.0; // user beta only on first update
            if(i_ == j_) {
               // diagonal block
               host_gemm(
                     OP_N, OP_T, u_ncol, u_ncol, cdata_[elim_col].nelim,
                     -1.0, &ld[ncol()], ldld,
                     &jsrc.aval_[ncol()], lda_,
                     beta, upd, ldupd
                     );
            } else {
               // off-diagonal block
               T* upd_ij =
                  &upd[(i_-calc_nblk(n_,block_size_))*block_size_+u_ncol];
               host_gemm(
                     OP_N, OP_T, nrow(), u_ncol, cdata_[elim_col].nelim,
                     -1.0, &ld[rfrom], ldld, &jsrc.aval_[ncol()], lda_,
                     beta, upd_ij, ldupd
                     );
            }
         }
      } else {
         // Update to left of elim column (UpdateT)
         int elim_col = jsrc.i_;
         if(cdata_[elim_col].nelim == 0) return; // nothing to do
         int rfrom = (i_ <= elim_col) ? cdata_[i_].nelim : 0;
         int cfrom = (j_ <= elim_col) ? cdata_[j_].nelim : 0;
         int ldld = align_lda<T>(block_size_);
         T* ld = work.get_ptr<T>(block_size_*ldld);
         // NB: we use ld[rfrom] below so alignment matches that of aval[rfrom]
         if(isrc.j_==elim_col) {
            calcLD<OP_N>(
                  nrow()-rfrom, cdata_[elim_col].nelim,
                  &isrc.aval_[rfrom], lda_,
                  cdata_[elim_col].d, &ld[rfrom], ldld
                  );
         } else {
            calcLD<OP_T>(
                  nrow()-rfrom, cdata_[elim_col].nelim, &
                  isrc.aval_[rfrom*lda_], lda_,
                  cdata_[elim_col].d, &ld[rfrom], ldld
                  );
         }
         host_gemm(
               OP_N, OP_N, nrow()-rfrom, ncol()-cfrom, cdata_[elim_col].nelim,
               -1.0, &ld[rfrom], ldld, &jsrc.aval_[cfrom*lda_], lda_,
               1.0, &aval_[cfrom*lda_+rfrom], lda_
               );
      }
   }

   void form_contrib(Block const& isrc, Block const& jsrc, Workspace& work, double beta, T* upd_ij, int ldupd) {
      int elim_col = isrc.j_;
      int ldld = align_lda<T>(block_size_);
      T* ld = work.get_ptr<T>(block_size_*ldld);
      calcLD<OP_N>(
            nrow(), cdata_[elim_col].nelim, isrc.aval_, lda_,
            cdata_[elim_col].d, ld, ldld
            );
      // User-supplied beta only on first update; otherwise 1.0
      T rbeta = (cdata_[elim_col].first_elim) ? beta : 1.0;
      int blkn = get_nrow(j_); // nrow not ncol as we're on contrib
      host_gemm(
            OP_N, OP_T, nrow(), blkn, cdata_[elim_col].nelim,
            -1.0, ld, ldld, jsrc.aval_, lda_,
            rbeta, upd_ij, ldupd
            );
   }

   /** Returns true if block contains NaNs (debug only). */
   bool isnan(int elim_col=-1) const {
      int m = (i_==elim_col) ? cdata_[i_].get_npass() : nrow();
      int n = (j_==elim_col) ? cdata_[j_].get_npass() : ncol();
      for(int j=0; j<n; ++j)
      for(int i=((i_==j_)?j:0); i<m; ++i) {
         if(std::isnan(aval_[j*lda_+i])) {
            printf("%d, %d is nan\n", i, j);
            return true;
         }
         if(!std::isfinite(aval_[j*lda_+i])) {
            printf("%d, %d is inf\n", i, j);
            return true;
         }
      }
      return false;
   }

   /** Prints block (debug only) */
   void print() const {
      printf("Block %d, %d (%d x %d):\n", i_, j_, nrow(), ncol());
      for(int i=0; i<nrow(); ++i) {
         printf("%d:", i);
         for(int j=0; j<ncol(); ++j)
            printf(" %e", aval_[j*lda_+i]);
         printf("\n");
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
   ColumnData<T,IntAlloc>& cdata_; //< global column data array
   T* aval_;
};

template<typename T,
         int BLOCK_SIZE,
         typename Backup,
         bool use_tasks,
         bool debug,
         typename Allocator
         >
class LDLT {
   typedef typename std::allocator_traits<Allocator>::template rebind_alloc<int> IntAlloc;
   typedef typename std::allocator_traits<Allocator>::template rebind_alloc<T> TAlloc;
private:
   /** Performs LDL^T factorization with block pivoting. Detects failure
    *  and aborts only column if an a posteori pivot test fails. */
   static
   int run_elim_pivoted(int const m, int const n, int* perm, T* a,
         int const lda, T* d, ColumnData<T,IntAlloc>& cdata, Backup& backup,
         struct cpu_factor_options const& options, int const block_size,
         T const beta, T* upd, int const ldupd, std::vector<Workspace>& work,
         Allocator const& alloc) {
      typedef Block<T, BLOCK_SIZE, IntAlloc> BlockSpec;
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
            shared(a, perm, backup, cdata, next_elim, d, \
                   options, work, alloc) \
            depend(inout: a[blk*block_size*lda+blk*block_size:1]) \
            depend(inout: perm[blk*block_size:1]) \
            if(use_tasks && mblk>1)
         {
#ifdef PROFILE
            Profile::Task task("TA_LDLT_DIAG", omp_get_thread_num());
#endif
            if(debug) printf("Factor(%d)\n", blk);
            BlockSpec dblk(blk, blk, m, n, cdata, a, lda, block_size);
            // Store a copy for recovery in case of a failed column
            dblk.backup(backup);
            // Perform actual factorization
            int nelim = dblk.template factor<Allocator>(
                  next_elim, perm, d, options, work, alloc
                  );
            // Init threshold check (non locking => task dependencies)
            cdata[blk].init_passed(nelim);
#ifdef PROFILE
            if(use_tasks) task.done();
#endif
         }
         
         // Loop over off-diagonal blocks applying pivot
         for(int jblk=0; jblk<blk; jblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, jblk) \
               shared(a, backup, cdata, options) \
               depend(in: a[blk*block_size*lda+blk*block_size:1]) \
               depend(inout: a[jblk*block_size*lda+blk*block_size:1]) \
               depend(in: perm[blk*block_size:1]) \
               if(use_tasks && mblk>1)
            {
#ifdef PROFILE
               Profile::Task task("TA_LDLT_APPLY", omp_get_thread_num());
#endif
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
#ifdef PROFILE
               if(use_tasks) task.done();
#endif
            }
         }
         for(int iblk=blk+1; iblk<mblk; iblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, iblk) \
               shared(a, backup, cdata, options) \
               depend(in: a[blk*block_size*lda+blk*block_size:1]) \
               depend(inout: a[blk*block_size*lda+iblk*block_size:1]) \
               depend(in: perm[blk*block_size:1]) \
               if(use_tasks && mblk>1)
            {
#ifdef PROFILE
               Profile::Task task("TA_LDLT_APPLY", omp_get_thread_num());
#endif
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
#ifdef PROFILE
               if(use_tasks) task.done();
#endif
            }
         }

         // Adjust column once all applys have finished and we know final
         // number of passed columns.
         #pragma omp task default(none) \
            firstprivate(blk) \
            shared(cdata, next_elim) \
            depend(inout: perm[blk*block_size:1]) \
            if(use_tasks && mblk>1)
         {
#ifdef PROFILE
            Profile::Task task("TA_LDLT_ADJUST", omp_get_thread_num());
#endif
            if(debug) printf("Adjust(%d)\n", blk);
            cdata[blk].adjust(next_elim);
#ifdef PROFILE
            if(use_tasks) task.done();
#endif
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
                  shared(a, cdata, backup, work)\
                  depend(inout: a[jblk*block_size*lda+iblk*block_size:1]) \
                  depend(in: perm[blk*block_size:1]) \
                  depend(in: a[jblk*block_size*lda+blk*block_size:1]) \
                  depend(in: a[adep_idx:1]) \
                  if(use_tasks && mblk>1)
               {
#ifdef PROFILE
                  Profile::Task task("TA_LDLT_UPDA", omp_get_thread_num());
#endif
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
                  ublk.update(isrc, jsrc, work[thread_num]);
#ifdef PROFILE
                  if(use_tasks) task.done();
#endif
               }
            }
         }
         for(int jblk=blk; jblk<nblk; jblk++) {
            for(int iblk=jblk; iblk<mblk; iblk++) {
               #pragma omp task default(none) \
                  firstprivate(blk, iblk, jblk) \
                  shared(a, cdata, backup, work, upd) \
                  depend(inout: a[jblk*block_size*lda+iblk*block_size:1]) \
                  depend(in: perm[blk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+iblk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+jblk*block_size:1]) \
                  if(use_tasks && mblk>1)
               {
#ifdef PROFILE
                  Profile::Task task("TA_LDLT_UPDA", omp_get_thread_num());
#endif
                  if(debug) printf("UpdateN(%d,%d,%d)\n", iblk, jblk, blk);
                  int thread_num = omp_get_thread_num();
                  BlockSpec ublk(iblk, jblk, m, n, cdata, a, lda, block_size);
                  BlockSpec isrc(iblk, blk, m, n, cdata, a, lda, block_size);
                  BlockSpec jsrc(jblk, blk, m, n, cdata, a, lda, block_size);
                  // If we're on the block col we've just eliminated, restore
                  // any failed cols and release resources storing backup
                  ublk.restore_if_required(backup, blk);
                  // Perform actual update
                  ublk.update(isrc, jsrc, work[thread_num],
                        beta, upd, ldupd);
#ifdef PROFILE
                  if(use_tasks) task.done();
#endif
               }
            }
         }

         // Handle update to contribution block, if required
         if(upd && mblk>nblk) {
            int uoffset = std::min(nblk*block_size, m) - n;
            T *upd2 = &upd[uoffset*(ldupd+1)];
            for(int jblk=nblk; jblk<mblk; ++jblk)
            for(int iblk=jblk; iblk<mblk; ++iblk) {
               T* upd_ij = &upd2[(jblk-nblk)*block_size*ldupd + 
                                 (iblk-nblk)*block_size];
               #pragma omp task default(none) \
                  firstprivate(iblk, jblk, blk, upd_ij) \
                  shared(a, upd2, cdata, work) \
                  depend(inout: upd_ij[0:1]) \
                  depend(in: perm[blk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+iblk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+jblk*block_size:1]) \
                  if(use_tasks && mblk>1)
               {
#ifdef PROFILE
                  Profile::Task task("TA_LDLT_UPDC", omp_get_thread_num());
#endif
                  if(debug) printf("FormContrib(%d,%d,%d)\n", iblk, jblk, blk);
                  int thread_num = omp_get_thread_num();
                  BlockSpec ublk(iblk, jblk, m, n, cdata, a, lda, block_size);
                  BlockSpec isrc(iblk, blk, m, n, cdata, a, lda, block_size);
                  BlockSpec jsrc(jblk, blk, m, n, cdata, a, lda, block_size);
                  ublk.form_contrib(
                        isrc, jsrc, work[thread_num], beta, upd_ij, ldupd
                        );
#ifdef PROFILE
                  if(use_tasks) task.done();
#endif
               }
            }
         }
      }
      if(use_tasks && mblk > 1) {
         // We only need a taskwait here if we've launched any subtasks...
         // NB: we don't use taskgroup as it doesn't support if()
         #pragma omp taskwait
      }

      /*if(debug) {
         printf("PostElim:\n");
         print_mat(mblk, nblk, m, n, blkdata, cdata, lda);
      }*/

      return next_elim;
   }

   /** Performs LDL^T factorization assuming everything works. Detects failure
    *  and aborts entire thing if a posteori pivot test fails. */
   static
   int run_elim_unpivoted(int const m, int const n, int* perm, T* a,
         int const lda, T* d, ColumnData<T,IntAlloc>& cdata, Backup& backup,
         struct cpu_factor_options const& options, int const block_size,
         T const beta, T* upd, int const ldupd, std::vector<Workspace>& work,
         Allocator const& alloc) {
      typedef Block<T, BLOCK_SIZE, IntAlloc> BlockSpec;
      //printf("ENTRY %d %d vis %d %d %d\n", m, n, mblk, nblk, block_size);

      int const nblk = calc_nblk(n, block_size);
      int const mblk = calc_nblk(m, block_size);

      /* Setup */
      int next_elim = 0;

      /* Inner loop - iterate over block columns */
      #pragma omp taskgroup
      for(int blk=0; blk<nblk; blk++) {
         /*if(debug) {
            printf("Bcol %d:\n", blk);
            print_mat(mblk, nblk, m, n, blkdata, cdata, lda);
         }*/

         // Factor diagonal
         #pragma omp task default(none) \
            firstprivate(blk) \
            shared(a, perm, backup, cdata, next_elim, d, \
                   options, work, alloc) \
            depend(inout: a[blk*block_size*lda+blk*block_size:1]) \
            if(use_tasks && mblk>1)
         {
#ifdef PROFILE
            Profile::Task task("TA_LDLT_DIAG", omp_get_thread_num());
#endif
            if(debug) printf("Factor(%d)\n", blk);
            BlockSpec dblk(blk, blk, m, n, cdata, a, lda, block_size);
            // On first access to this block, store copy in case of failure
            if(blk==0) dblk.backup(backup);
            // Perform actual factorization
            int nelim = dblk.template factor<Allocator>(
                  next_elim, perm, d, options, work, alloc
                  );
            if(nelim < get_ncol(blk, n, block_size)) {
               #pragma omp cancel taskgroup
            }
            cdata[blk].first_elim = (blk==0);
            cdata[blk].init_passed(1); // diagonal block has passed
            next_elim += nelim; // we're assuming everything works
#ifdef PROFILE
            if(use_tasks) task.done();
#endif
         }
         
         // Loop over off-diagonal blocks applying pivot
         for(int jblk=0; jblk<blk; jblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, jblk) \
               shared(a, backup, cdata, options, work) \
               depend(in: a[blk*block_size*lda+blk*block_size:1]) \
               depend(inout: a[jblk*block_size*lda+blk*block_size:1]) \
               if(use_tasks && mblk>1)
            {
#ifdef PROFILE
               Profile::Task task("TA_LDLT_APPLY", omp_get_thread_num());
#endif
               if(debug) printf("ApplyT(%d,%d)\n", blk, jblk);
               int thread_num = omp_get_thread_num();
               BlockSpec dblk(blk, blk, m, n, cdata, a, lda, block_size);
               BlockSpec cblk(blk, jblk, m, n, cdata, a, lda, block_size);
               // Apply row permutation from factorization of dblk
               cblk.apply_rperm(work[thread_num]);
               // NB: no actual application of pivot must be done, as we are
               // assuming everything has passed...
#ifdef PROFILE
               if(use_tasks) task.done();
#endif
            }
         }
         for(int iblk=blk+1; iblk<mblk; iblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, iblk) \
               shared(a, backup, cdata, options, work) \
               depend(in: a[blk*block_size*lda+blk*block_size:1]) \
               depend(inout: a[blk*block_size*lda+iblk*block_size:1]) \
               if(use_tasks && mblk>1)
            {
#ifdef PROFILE
               Profile::Task task("TA_LDLT_APPLY", omp_get_thread_num());
#endif
               if(debug) printf("ApplyN(%d,%d)\n", iblk, blk);
               int thread_num = omp_get_thread_num();
               BlockSpec dblk(blk, blk, m, n, cdata, a, lda, block_size);
               BlockSpec rblk(iblk, blk, m, n, cdata, a, lda, block_size);
               // On first access to this block, store copy in case of failure
               if(blk==0) rblk.backup(backup);
               // Apply column permutation from factorization of dblk
               rblk.apply_cperm(work[thread_num]);
               // Perform elimination and determine number of rows in block
               // passing a posteori threshold pivot test
               int blkpass = rblk.apply_pivot_app(dblk, options.u, options.small);
               // Update column's passed pivot count
               if(cdata[blk].test_fail(blkpass)) {
                  #pragma omp cancel taskgroup
               }
#ifdef PROFILE
               if(use_tasks) task.done();
#endif
            }
         }

         // Update uneliminated columns
         // Column blk only needed if upd is present
         int jsa = (upd) ? blk : blk + 1;
         for(int jblk=jsa; jblk<nblk; jblk++) {
            for(int iblk=jblk; iblk<mblk; iblk++) {
               #pragma omp task default(none) \
                  firstprivate(blk, iblk, jblk) \
                  shared(a, cdata, backup, work, upd) \
                  depend(inout: a[jblk*block_size*lda+iblk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+iblk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+jblk*block_size:1]) \
                  if(use_tasks && mblk>1)
               {
#ifdef PROFILE
                  Profile::Task task("TA_LDLT_UPDA", omp_get_thread_num());
#endif
                  if(debug) printf("UpdateN(%d,%d,%d)\n", iblk, jblk, blk);
                  int thread_num = omp_get_thread_num();
                  BlockSpec ublk(iblk, jblk, m, n, cdata, a, lda, block_size);
                  BlockSpec isrc(iblk, blk, m, n, cdata, a, lda, block_size);
                  BlockSpec jsrc(jblk, blk, m, n, cdata, a, lda, block_size);
                  // On first access to this block, store copy in case of fail
                  if(blk==0 && jblk!=blk) ublk.backup(backup);
                  // Actual update
                  ublk.update(isrc, jsrc, work[thread_num], beta, upd, ldupd);
#ifdef PROFILE
                  if(use_tasks) task.done();
#endif
               }
            }
         }

         // Handle update to contribution block, if required
         if(upd && mblk>nblk) {
            int uoffset = std::min(nblk*block_size, m) - n;
            T *upd2 = &upd[uoffset*(ldupd+1)];
            for(int jblk=nblk; jblk<mblk; ++jblk)
            for(int iblk=jblk; iblk<mblk; ++iblk) {
               T* upd_ij = &upd2[(jblk-nblk)*block_size*ldupd + 
                                 (iblk-nblk)*block_size];
               #pragma omp task default(none) \
                  firstprivate(iblk, jblk, blk, upd_ij) \
                  shared(a, upd2, cdata, work) \
                  depend(inout: upd_ij[0:1]) \
                  depend(in: a[blk*block_size*lda+iblk*block_size:1]) \
                  depend(in: a[blk*block_size*lda+jblk*block_size:1]) \
                  if(use_tasks && mblk>1)
               {
#ifdef PROFILE
                  Profile::Task task("TA_LDLT_UPDC", omp_get_thread_num());
#endif
                  if(debug) printf("FormContrib(%d,%d,%d)\n", iblk, jblk, blk);
                  int thread_num = omp_get_thread_num();
                  BlockSpec ublk(iblk, jblk, m, n, cdata, a, lda, block_size);
                  BlockSpec isrc(iblk, blk, m, n, cdata, a, lda, block_size);
                  BlockSpec jsrc(jblk, blk, m, n, cdata, a, lda, block_size);
                  ublk.form_contrib(
                        isrc, jsrc, work[thread_num], beta, upd_ij, ldupd
                        );
#ifdef PROFILE
                  if(use_tasks) task.done();
#endif
               }
            }
         }
      }
      // FIXME: ...
      /*if(use_tasks && mblk > 1) {
         // We only need a taskwait here if we've launched any subtasks...
         // NB: we don't use taskgroup as it doesn't support if()
         #pragma omp taskwait
      }*/

      /*if(debug) {
         printf("PostElim:\n");
         print_mat(mblk, nblk, m, n, blkdata, cdata, lda);
      }*/

      return cdata.calc_nelim(m);
   }

   /** Restore matrix to original state prior to aborted factorization */
   static
   void restore(int const m, int const n, int* perm, T* a, int const lda, T* d,
         ColumnData<T,IntAlloc>& cdata, Backup& backup, int const* old_perm,
         int const block_size) {
      typedef Block<T, BLOCK_SIZE, IntAlloc> BlockSpec;

      int const nblk = calc_nblk(n, block_size);
      int const mblk = calc_nblk(m, block_size);

      /* Restore perm */
      for(int i=0; i<n; ++i)
         perm[i] = old_perm[i];

      /* Restore a */
      for(int jblk=0; jblk<nblk; ++jblk) {
         for(int iblk=jblk; iblk<mblk; ++iblk) {
            #pragma omp task default(none) \
               firstprivate(iblk, jblk) \
               shared(a, backup, cdata)
            {
               BlockSpec rblk(iblk, jblk, m, n, cdata, a, lda, block_size);
               rblk.full_restore(backup);
            }
         }
      }
      #pragma omp taskwait
   }


   static
   void print_mat(int m, int n, const int *perm, std::vector<bool> const& eliminated, const T *a, int lda) {
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
   int factor(int m, int n, int *perm, T *a, int lda, T *d, Backup& backup, struct cpu_factor_options const& options, PivotMethod pivot_method, int block_size, T beta, T* upd, int ldupd, std::vector<Workspace>& work, Allocator const& alloc=Allocator()) {
      /* Sanity check arguments */
      if(m < n) return -1;
      if(lda < n) return -4;

      /* Initialize useful quantities: */
      int nblk = calc_nblk(n, block_size);
      int mblk = calc_nblk(m, block_size);

      /* Temporary workspaces */
      ColumnData<T, IntAlloc> cdata(n, block_size, IntAlloc(alloc));
#ifdef PROFILE
      Profile::setNullState(omp_get_thread_num());
#endif

      /* Main loop
       *    - Each pass leaves any failed pivots in place and keeps everything
       *      up-to-date.
       *    - If no pivots selected across matrix, perform swaps to get large
       *      entries into diagonal blocks
       */
      int num_elim;
      if(pivot_method == PivotMethod::app_aggressive) {
         // Take a copy of perm
         typedef std::allocator_traits<IntAlloc> IATraits;
         IntAlloc intAlloc(alloc);
         int* perm_copy = IATraits::allocate(intAlloc, n);
         for(int i=0; i<n; ++i) {
            perm_copy[i] = perm[i];
         }
         if(beta!=0.0) {
            throw std::runtime_error(
                  "run_elim_unpivoted currently only supports beta=0.0"
                  );
         }
         // Run the elimination
         num_elim = run_elim_unpivoted(
               m, n, perm, a, lda, d, cdata, backup, options,
               block_size, beta, upd, ldupd, work, alloc
               );
         if(num_elim < n) {
#ifdef PROFILE
            {
               char buffer[200];
               snprintf(buffer, 200, "tpp-aggressive failed at %d / %d\n",
                        num_elim, n);
               Profile::addEvent("EV_AGG_FAIL", omp_get_thread_num(), buffer);
            }
#endif
            // Factorization ecountered a pivoting failure.
            // Rollback to known good state
            restore(
                  m, n, perm, a, lda, d, cdata, backup, perm_copy, block_size
                  );
            // Factorize more carefully
            num_elim =
               LDLT<T, BLOCK_SIZE, Backup, use_tasks, debug, Allocator>
               ::factor(
                     m, n, perm, a, lda, d, backup, options,
                     PivotMethod::app_block, block_size, beta,
                     upd, ldupd, work, alloc
                     );
         }
         IATraits::deallocate(intAlloc, perm_copy, n);
         return num_elim;
      } else {
         num_elim = run_elim_pivoted(
               m, n, perm, a, lda, d, cdata, backup, options,
               block_size, beta, upd, ldupd, work, alloc
               );
         backup.release_all_memory(); // we're done with it now, but we want
                                      // the memory back for reuse before we
                                      // get it automatically when it goes out
                                      // of scope.

         // Permute failed entries to end
#ifdef PROFILE
         Profile::Task task_post("TA_LDLT_POST", omp_get_thread_num());
#endif
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
#ifdef PROFILE
         task_post.done();
#endif
      }

      if(debug) {
         std::vector<bool> eliminated(n);
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
size_t ldlt_app_factor_mem_required(int m, int n, int block_size) {
   int const align = 32;
   return align_lda<T>(m) * n * sizeof(T) + align; // CopyBackup
}

template<typename T, typename Allocator>
int ldlt_app_factor(int m, int n, int* perm, T* a, int lda, T* d, T beta, T* upd, int ldupd, struct cpu_factor_options const& options, std::vector<Workspace>& work, Allocator const& alloc) {
   // If we've got a tall and narrow node, adjust block size so each block
   // has roughly blksz**2 entries
   // FIXME: Decide if this reshape is actually useful, given it will generate
   //        a lot more update tasks instead?
   int outer_block_size = options.cpu_task_block_size;
   /*if(n < outer_block_size) {
       outer_block_size = int((long(outer_block_size)*outer_block_size) / n);
   }*/

#ifdef PROFILE
   Profile::setState("TA_MISC1", omp_get_thread_num());
#endif

   // Template parameters and workspaces
   bool const debug = false;
   //PoolBackup<T, Allocator> backup(m, n, outer_block_size, alloc);
   CopyBackup<T, Allocator> backup(m, n, outer_block_size, alloc);

   // Actual call
   bool const use_tasks = true;
   return LDLT
      <T, INNER_BLOCK_SIZE, CopyBackup<T,Allocator>, use_tasks, debug,
       Allocator>
      ::factor(
            m, n, perm, a, lda, d, backup, options, options.pivot_method,
            outer_block_size, beta, upd, ldupd, work, alloc
            );
}
template int ldlt_app_factor<double, BuddyAllocator<double,std::allocator<double>>>(int, int, int*, double*, int, double*, double, double*, int, struct cpu_factor_options const&, std::vector<Workspace>&, BuddyAllocator<double,std::allocator<double>> const& alloc);

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
