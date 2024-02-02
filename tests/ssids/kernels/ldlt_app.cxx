/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * Licence: BSD licence, see LICENCE file for details
 *
 */
#include "ldlt_app.hxx"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>

#include "AlignedAllocator.hxx"
#include "framework.hxx"
#include "ssids/cpu/kernels/wrappers.hxx"
#include "ssids/cpu/kernels/ldlt_app.cxx" // .cxx as we need internal namespace
#include "ssids/cpu/kernels/ldlt_tpp.hxx"
#include "ssids/cpu/cpu_iface.hxx"

using namespace spral::ssids::cpu;

template<typename T>
void print_d(int n, T *d) {
   bool first = true;
   for(int i=0; i<n; i++) {
      if(d[2*i]==std::numeric_limits<T>::infinity() && d[2*i+1]!=0.0) {
         // Second half of 2x2 pivot: don't print brackets
         std::cout << " ";
      } else if(first) {
         std::cout << "(";
      } else {
         std::cout << ") (";
      }
      if(d[2*i] == std::numeric_limits<T>::infinity() && d[2*i+1]!=0.0) std::cout << std::setw(8) << " ";
      else              std::cout << std::setw(8) << d[2*i];
      first = false;
   }
   std::cout << ")" << std::endl;
   first = true;
   for(int i=0; i<n; i++) {
      if(d[2*i]==std::numeric_limits<T>::infinity() && d[2*i+1]!=0.0) {
         // Second half of 2x2 pivot: don't print brackets
         std::cout << " ";
      } else if(first) {
         std::cout << "(";
      } else {
         std::cout << ") (";
      }
      if(d[2*i+1] == 0.0) std::cout << std::setw(8) << " ";
      else                std::cout << std::setw(8) << d[2*i+1];
      first = false;
   }
   std::cout << ")" << std::endl;
}

template<typename T>
void solve(int m, int n, const int *perm, const T *l, int ldl, const T *d, const T *b, T *x) {
   for(int i=0; i<m; i++) x[i] = b[perm[i]];
   // Fwd slv
   ldlt_app_solve_fwd(m, n, l, ldl, 1, x, m);
   ldlt_app_solve_fwd(m-n, m-n, &l[n*(ldl+1)], ldl, 1, &x[n], m);
   // Diag slv
   ldlt_app_solve_diag(n, d, 1, x, m);
   ldlt_app_solve_diag(m-n, &d[2*n], 1, &x[n], m);
   // Bwd slv
   ldlt_app_solve_bwd(m-n, m-n, &l[n*(ldl+1)], ldl, 1, &x[n], m);
   ldlt_app_solve_bwd(m, n, l, ldl, 1, x, m);
   // Undo permutation
   T *temp = new T[m];
   for(int i=0; i<m; i++) temp[i] = x[i];
   for(int i=0; i<m; i++) x[perm[i]] = temp[i];
   // Free mem
   delete[] temp;
}

/// Makes a (symmetric, half storage) matrix singular by making col2 an
/// appropriate multiple of col1
template <typename T>
void make_singular(int n, int col1, int col2, T *a, int lda) {
   T *col = new T[n];

   T a11 = a[col1*(lda+1)];
   T a21 = (col1 < col2) ? a[col1*lda + col2]
                         : a[col2*lda + col1];
   T scal = a21 / a11;

   // Read col1 and double it
   for(int i=0; i<col1; i++)
      col[i] = scal*a[i*lda+col1];
   for(int i=col1; i<n; i++)
      col[i] = scal*a[col1*lda+i];

   // Store col to col2
   for(int i=0; i<col2; i++)
      a[i*lda+col2] = col[i];
   for(int i=col2; i<n; i++)
      a[col2*lda+i] = col[i];

   // Free mem
   delete[] col;
}

/// Makes a specified diagonal block of a matrix singular by making first and
//  last columns linearlly dependent
template <typename T, int BLOCK_SIZE>
void make_dblk_singular(int blk, int nblk, T *a, int lda) {
   int col1 = 0;
   int col2 = BLOCK_SIZE-1;
   T *adiag = &a[blk*BLOCK_SIZE*(lda+1)];
   make_singular(BLOCK_SIZE, col1, col2, adiag, lda);
}

// Pick n/8 random rows and multiply by 1000. Then do the same for n/8 random entries.
template <typename T, int BLOCK_SIZE>
void cause_delays(int n, T *a, int lda) {
   int nsing = n/8;
   if(nsing==0) nsing=1;
   for(int i=0; i<nsing; i++) {
      // Add a row of oversized values
      int idx = n*((float) rand())/RAND_MAX;
      if(i==0 && n>BLOCK_SIZE && idx<BLOCK_SIZE) idx += BLOCK_SIZE;
      idx = std::min(idx, n);
      for(int c=0; c<idx; c++)
         a[c*lda+idx] *= 1000;
      for(int r=idx; r<n; r++)
         a[idx*lda+r] *= 1000;
      int row = n*((float) rand())/RAND_MAX;
      int col = n*((float) rand())/RAND_MAX;
      if(row > col) a[col*lda+row] *= 1000;
      else          a[row*lda+col] *= 1000;
   }
}

// Reorders a into supplied permuation such that row/col i of the new a
// contains row/col perm[i] of the original a
template <typename T>
void reorder(int n, const int *perm, T *a, int lda) {
   T *a2 = new T[n*n];
   // Copy to a2 with permutation
   for(int j=0; j<n; j++) {
      int c = perm[j];
      for(int i=j; i<n; i++) {
         int r = perm[i];
         a2[j*n+i] = (r<c) ? a[r*lda+c]
                           : a[c*lda+r];
      }
   }
   // Copy back to a
   for(int j=0; j<n; j++)
   for(int i=j; i<n; i++)
      a[j*n+i] = a2[j*n+i];
   delete[] a2;
}

/// Update A22 -= A_21 D_11 A_21^T
template <typename T>
void do_update(int n, int k, T *a22, const T *a21, int lda, const T* d) {
   // Form A_21 D_11
   T *ad21 = new T[n*k];
   for(int j=0; j<k;) {
      if(j+1==k || std::isfinite(d[2*j+2])) {
         // 1x1 pivot
         // (Actually stored as D^-1 so need to invert it again)
         if(d[2*j] == 0.0) {
            // Handle zero pivots with care
            for(int i=0; i<n; i++) {
               ad21[j*n+i] = 0.0;
            }
         } else {
            // Standard 1x1 pivot
            T d11 = 1/d[2*j];
            // And calulate ad21
            for(int i=0; i<n; i++) {
               ad21[j*n+i] = d11*a21[j*lda+i];
            }
         }
         // Increment j
         j++;
      } else {
         // 2x2 pivot
         // (Actually stored as D^-1 so need to invert it again)
         T di11 = d[2*j]; T di21 = d[2*j+1]; T di22 = d[2*j+3];
         T det = di11*di22 - di21*di21;
         T d11 = di22 / det; T d21 = -di21 / det; T d22 = di11 / det;
         // And calulate ad21
         for(int i=0; i<n; i++) {
            ad21[j*n+i]     = d11*a21[j*lda+i] + d21*a21[(j+1)*lda+i];
            ad21[(j+1)*n+i] = d21*a21[j*lda+i] + d22*a21[(j+1)*lda+i];
         }
         // Increment j
         j += 2;
      }
   }

   /*printf("a21:\n");
   for(int i=0; i<n; i++) {
      for(int j=0; j<k; j++)
         printf(" %le", a21[j*lda+i]);
      printf("\n");
   }
   printf("ad21:\n");
   for(int i=0; i<n; i++) {
      for(int j=0; j<k; j++)
         printf(" %le", ad21[j*n+i]);
      printf("\n");
   }
   printf("a22:\n");
   for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++)
         printf(" %le", a22[j*lda+i]);
      printf("\n");
   }*/

   // Perform actual update
   host_gemm<T>(OP_N, OP_T, n, n, k, -1.0, ad21, n, a21, lda, 1.0, a22, lda);

   // Free memory
   delete[] ad21;
}

/// Permutes rows of a as per perm such that row i becomes row (perm[i]-offset)
template<typename T>
void permute_rows(int n, int k, const int *reord, int *perm, T *a, int lda) {
   // Copy a and perm without permutation
   T *a2 = new T[n*k];
   for(int j=0; j<k; j++)
   for(int i=0; i<n; i++)
      a2[j*n+i] = a[j*lda+i];
   int *perm2 = new int[n];
   for(int i=0; i<n; i++)
      perm2[i] = perm[i];

   // Copy back with permutation
   for(int j=0; j<k; j++)
   for(int i=0; i<n; i++) {
      int row = reord[i];
      a[j*lda+i] = a2[j*n+row];
   }
   for(int i=0; i<n; i++) {
      int row = reord[i];
      perm[i] = perm2[row];
   }

   // Free memory
   delete[] a2;
   delete[] perm2;
}

template <typename T,
          int BLOCK_SIZE
          >
void modify_test_matrix(bool singular, bool delays, bool dblk_singular, int m, int n, T *a, int lda) {
   int mblk = m / BLOCK_SIZE;
   int nblk = n / BLOCK_SIZE;
   if(delays)
      cause_delays<T,BLOCK_SIZE>(m, a, lda);
   if(dblk_singular) {
      int blk = nblk * ((float) rand())/RAND_MAX;
      make_dblk_singular<T, BLOCK_SIZE>(blk, mblk, a, lda);
   }
   if(n>1 && singular) {
      int col1 = n * ((float) rand())/RAND_MAX;
      int col2 = col1;
      while(col1 == col2)
         col2 = n * ((float) rand())/RAND_MAX;
      make_singular<T>(m, col1, col2, a, lda);
   }
}

template <typename T,
          int INNER_BLOCK_SIZE,
          bool aggressive, // Use Cholesky-like app pattern
          bool debug // Switch on debugging output
          >
int ldlt_test(T u, T small, bool delays, bool singular, bool dblk_singular, int m, int n, int outer_block_size=INNER_BLOCK_SIZE, int test=0, int seed=0) {
   // Note: We generate an m x m test matrix, then factor it as an
   // m x n matrix followed by an (m-n) x (m-n) matrix [give or take delays]
   bool failed = false;

   // Generate test matrix
   int lda = align_lda<T>(m);
   double* a = new T[m*lda];
   gen_sym_indef(m, a, lda);
   modify_test_matrix<T, INNER_BLOCK_SIZE>(
         singular, delays, dblk_singular, m, n, a, lda
         );

   // Generate a RHS based on x=1, b=Ax
   T *b = new T[m];
   gen_rhs(m, a, lda, b);

   // Print out matrices if requested
   if(debug) {
      std::cout << "A:" << std::endl;
      print_mat("%10.2e", m, a, lda);
   }

   // Setup options
   struct cpu_factor_options options;
   options.action = true;
   options.multiplier = 2.0;
   options.small = small;
   options.u = u;
   options.print_level = 0;
   options.small_subtree_threshold = 100*100*100;
   options.cpu_block_size = 256;
   options.pivot_method = (aggressive) ? PivotMethod::app_aggressive
                                       : PivotMethod::app_block;

   // Factorize using main routine
   spral::test::AlignedAllocator<T> allocT;
   T *l = allocT.allocate(m*lda);
   memcpy(l, a, m*lda*sizeof(T)); // Copy a to l
   int *perm = new int[m];
   for(int i=0; i<m; i++) perm[i] = i;
   T *d = new T[2*m];
   // First m x n matrix
   CopyBackup<T> backup(m, n, outer_block_size);
   std::vector<Workspace> work;
   const int SSIDS_PAGE_SIZE = 8*1024*1024; // 8 MB
   for(int i=0; i<omp_get_num_threads(); ++i)
      work.emplace_back(SSIDS_PAGE_SIZE);
   int const use_tasks = true;
   int q1 = LDLT
         <T, INNER_BLOCK_SIZE, CopyBackup<T>, use_tasks, debug>
         ::factor(
            m, n, perm, l, lda, d, backup, options, options.pivot_method,
            outer_block_size, 0.0, nullptr, 0, work
            );
   if(debug) {
      std::cout << "FIRST FACTOR CALL ELIMINATED " << q1 << " of " << n << " pivots" << std::endl;
      std::cout << "L after first elim:" << std::endl;
      print_mat("%10.2e", m, l, lda, perm);
      std::cout << "D:" << std::endl;
      print_d<T>(q1, d);
   }
   int q2 = 0;
   if(q1 < n) {
      // Finish off with simplistic kernel
      T *ld = new T[2*m];
      q1 += ldlt_tpp_factor(m-q1, n-q1, &perm[q1], &l[(q1)*(lda+1)], lda,
            &d[2*(q1)], ld, m, options.action, u, small, q1, &l[q1], lda);
      delete[] ld;
   }
   if(m > n) {
      // Apply outer product update
      do_update<T>(m-n, q1, &l[n*(lda+1)], &l[n], lda, d);
      // Second (m-n) x (m-n) matrix [but add delays if any]
      int *perm2 = new int[m-q1];
      for(int i=0; i<m-q1; i++)
         perm2[i] = i;
      CopyBackup<T> backup(m-q1, m-q1, outer_block_size);
      q2 = LDLT
         <T, INNER_BLOCK_SIZE, CopyBackup<T>, use_tasks, debug>
         ::factor(
            m-q1, m-q1, perm2, &l[q1*(lda+1)], lda, &d[2*q1], backup, options,
            options.pivot_method, outer_block_size, 0.0, nullptr, 0, work
            );
      // Permute rows of A_21 as per perm
      permute_rows(m-q1, q1, perm2, &perm[q1], &l[q1], lda);
      delete[] perm2;
      if(q1+q2 < m) {
         // Finish off with simplistic kernel
         T *ld = new T[2*m];
         q2 += ldlt_tpp_factor(m-q1-q2, m-q1-q2, &perm[q1+q2],
               &l[(q1+q2)*(lda+1)], lda, &d[2*(q1+q2)], ld, m, options.action,
               u, small, q1+q2, &l[q1+q2], lda);
         delete[] ld;
      }
   }
   EXPECT_EQ(m, q1+q2) << "(test " << test << " seed " << seed << ")" << std::endl;

   // Print out matrices if requested
   if(debug) {
      std::cout << "q1=" << q1 << " q2=" << q2 << std::endl;
      std::cout << "L:" << std::endl;
      print_mat("%10.2e", m, l, lda, perm);
      std::cout << "D:" << std::endl;
      print_d<T>(m, d);
   }

   // Perform solve
   T *soln = new T[m];
   solve(m, q1, perm, l, lda, d, b, soln);
   if(debug) {
      printf("soln = ");
      for(int i=0; i<m; i++) printf(" %le", soln[i]);
      printf("\n");
   }

   // Check residual
   T bwderr = backward_error(m, a, lda, b, 1, soln, m);
   if(debug) printf("bwderr = %le\n", bwderr);
   EXPECT_LE(bwderr, 1e-13) << "(test " << test << " seed " << seed << ")" << std::endl;

   // Cleanup memory
   delete[] a; allocT.deallocate(l, m*lda);
   delete[] b;
   delete[] perm;
   delete[] d; delete[] soln;

   return (failed) ? -1 : 0;
}

template <typename T>
void print_mat (int n, int *perm, T *a, int lda) {
   for(int i=0; i<n; i++) {
      printf("%d:", perm[i]);
      for(int j=0; j<n; j++)
         printf(" %le", a[j*lda+i]);
      printf("\n");
   }
}

template <typename T,
          int BLOCK_SIZE,
          bool aggressive, // Use Cholesky-like app pattern
          int ntest, // Number of tests
          bool debug // Switch on debugging output
          >
int ldlt_torture_test(T u, T small, int m, int n) {
   for(int test=0; test<ntest; test++) {
      // Record seed we're using
      unsigned int seed = rand();
      //seed = 1872440417;
      srand(seed);

      // Matrix has:
      //    70% chance of getting delays
      //    20% chance of being made singular
      //    10% chance of a singular diagonal block
      bool delays = ( ((float) rand())/RAND_MAX < 0.7 );
      bool singular = ( ((float) rand())/RAND_MAX < 0.2 );
      bool dblk_singular = ( ((float) rand())/RAND_MAX < 0.1 );

      // Output debug info
      if(debug) {
         std::cout << "##########################################" << std::endl;
         std::cout << "# Test " << test << " singular:" << (singular?"T":"F") << " delays:" << (delays?"T":"F") << std::endl;
         std::cout << "# Random seed " << seed << std::endl;
         std::cout << "##########################################" << std::endl;
      }

      int err = ldlt_test<T, BLOCK_SIZE, aggressive, debug>(u, small, delays, singular, dblk_singular, m, n, BLOCK_SIZE, test, seed);
      if(err!=0) return err;
   }

   return 0; // Success
}

int run_ldlt_app_tests() {
   int nerr = 0;

   /* Simple tests, single level blocking, rectangular, app_block */
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 2, 1)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 4, 2)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 5, 3)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 8, 2)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, false, 64, 24)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, false, 23, 9)
      ));

   /* Simple tests, single level blocking, rectangular, app_aggressive */
   TEST((
      ldlt_test<double, 2, true, false> (0.01, 1e-20, false, false, false, 2, 1)
      ));
   TEST((
      ldlt_test<double, 2, true, false> (0.01, 1e-20, false, false, false, 4, 2)
      ));
   TEST((
      ldlt_test<double, 2, true, false> (0.01, 1e-20, false, false, false, 5, 3)
      ));
   TEST((
      ldlt_test<double, 2, true, false> (0.01, 1e-20, false, false, false, 8, 2)
      ));
   TEST((
      ldlt_test<double, 8, true, false> (0.01, 1e-20, false, false, false, 64, 24)
      ));
   TEST((
      ldlt_test<double, 8, true, false> (0.01, 1e-20, false, false, false, 23, 9)
      ));

   /* Simple tests, two-level blocking, rectangular */
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 2, 1, 4)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 4, 2, 4)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 5, 3, 4)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 8, 2, 4)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, false, 23, 9, 32)
      ));
   TEST((
      ldlt_test<double, 4, false, false> (0.01, 1e-20, false, false, false, 32, 12, 8)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, false, 64, 24, 32)
      ));

   /* Simple tests, single-level blocking, rectangular, non-singular,
    * with delays */
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, true, false, false, 4, 2)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, true, false, false, 8, 2)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, true, false, false, 64, 24)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, true, false, false, 23, 9)
      ));

   /* Simple tests, single-level blocking, rectangular, non-singular,
    * with delays, app-aggressive */
   TEST((
      ldlt_test<double, 2, true, false> (0.01, 1e-20, true, false, false, 4, 2)
      ));
   TEST((
      ldlt_test<double, 2, true, false> (0.01, 1e-20, true, false, false, 8, 2)
      ));
   TEST((
      ldlt_test<double, 8, true, false> (0.01, 1e-20, true, false, false, 64, 24)
      ));
   TEST((
      ldlt_test<double, 8, true, false> (0.01, 1e-20, true, false, false, 23, 9)
      ));

   /* Simple tests, single-level blocking, square, non-singular, no delays */
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 1*2, 1*2)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 2*2, 2*2)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, false, false, false, 8*2, 8*2)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, false, 8*8, 8*8)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, false, 27, 27)
      ));

   /* Simple tests, single-level blocking, square, non-singular, with delays */
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, true, false, false, 1*2, 1*2)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, true, false, false, 2*2, 2*2)
      ));
   TEST((
      ldlt_test<double, 2, false, false> (0.01, 1e-20, true, false, false, 4*2, 4*2)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, true, false, false, 8*8, 8*8)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, true, false, false, 29, 29)
      ));

   /* Test edge case of singular diagonal blocks, square non-singular matrix */
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, true, 8*8, 8*8)
      ));
   TEST((
      ldlt_test<double, 8, false, false> (0.01, 1e-20, false, false, true, 33, 33)
      ));

   /* Torture tests */
   TEST((
      ldlt_torture_test<double, 16, false, 500, false> (0.01, 1e-20, 8*16, 8*16)
      ));
   TEST((
      ldlt_torture_test<double, 16, false, 500, false> (0.01, 1e-20, 8*16, 3*16)
      ));

   return nerr;
}
