/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * Licence: BSD licence, see LICENCE file for details
 *
 */
#include "ldlt_tpp.hxx"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>

#include "framework.hxx"
#include "ssids/cpu/kernels/wrappers.hxx"
#include "ssids/cpu/kernels/ldlt_tpp.hxx"

using namespace spral::ssids::cpu;

namespace {

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
   ldlt_tpp_solve_fwd(m, n, l, ldl, 1, x, m);
   ldlt_tpp_solve_fwd(m-n, m-n, &l[n*(ldl+1)], ldl, 1, &x[n], m);
   // Diag slv
   ldlt_tpp_solve_diag(n, d, x);
   ldlt_tpp_solve_diag(m-n, &d[2*n], &x[n]);
   // Bwd slv
   ldlt_tpp_solve_bwd(m-n, m-n, &l[n*(ldl+1)], ldl, 1, &x[n], m);
   ldlt_tpp_solve_bwd(m, n, l, ldl, 1, x, m);
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

// Pick n/8 random rows and multiply by 1000. Then do the same for n/8 random entries.
void cause_delays(int n, double *a, int lda) {
   int nsing = n/8;
   if(nsing==0) nsing=1;
   for(int i=0; i<nsing; i++) {
      // Add a row of oversized values
      int idx = n*((float) rand())/RAND_MAX;
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
      if(j+1<k && std::isinf(d[2*j+2])) {
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
      } else {
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

void modify_test_matrix(bool singular, bool delays, int m, int n, double *a, int lda) {
   if(delays)
      cause_delays(m, a, lda);
   if(singular && n!=1) {
      int col1 = n * ((float) rand())/RAND_MAX;
      int col2 = col1;
      while(col1 == col2)
         col2 = n * ((float) rand())/RAND_MAX;
      make_singular(m, col1, col2, a, lda);
   }
}

double find_l_abs_max(int n, double *a, int lda) {
   double best = 0.0;
   for(int c=0; c<n; ++c)
   for(int r=c; r<n; ++r)
      best = std::max(best, fabs(a[c*lda+r]));
   return best;
}

int ldlt_tpp_test(double u, double small, bool delays, bool singular, int m, int n, bool debug=false, int test=0, int seed=0) {
   bool failed = false;
   bool action = true; // Don't abort on singular matrices
   // Note: We generate an m x m test matrix, then factor it as an
   // m x n matrix followed by an (m-n) x (m-n) matrix [give or take delays]

   // Generate test matrix
   int lda = m;
   double* a = new double[m*lda];
   gen_sym_indef(m, a, lda);
   modify_test_matrix(singular, delays, m, n, a, lda);

   // Generate a RHS based on x=1, b=Ax
   double *b = new double[m];
   gen_rhs(m, a, lda, b);

   // Print out matrices if requested
   if(debug) {
      std::cout << "A:" << std::endl;
      print_mat("%10.2e", m, a, lda);
   }

   // Factorize using main routine
   double *l = new double[m*lda];
   memcpy(l, a, m*lda*sizeof(double)); // Copy a to l
   int *perm = new int[m];
   for(int i=0; i<m; i++) perm[i] = i;
   double *d = new double[2*m];
   double *work = new double[2*m];
   // First m x n matrix
   int q1 = ldlt_tpp_factor(m, n, perm, l, lda, d, work, m, action, u, small);
   if(debug) std::cout << "FIRST FACTOR CALL ELIMINATED " << q1 << " of " << n << " pivots" << std::endl;
   int q2 = 0;
   if(m > n) {
      // Apply outer product update
      do_update<double>(m-n, q1, &l[n*(lda+1)], &l[n], lda, d);
      // Second (m-n) x (m-n) matrix [but add delays if any]
      q2 = ldlt_tpp_factor(m-q1, m-q1, &perm[q1], &l[q1*(lda+1)], lda, &d[2*q1], work, m, action, u, small, q1, &l[q1], lda);
   }
   EXPECT_EQ(m, q1+q2) << "(test " << test << " seed " << seed << ")" << std::endl;
   EXPECT_LE(find_l_abs_max(m, l, lda), 1.0/u) << "(test " << test << " seed " << seed << ")" << std::endl;

   // Print out matrices if requested
   if(debug) {
      std::cout << "q1=" << q1 << " q2=" << q2 << std::endl;
      std::cout << "L:" << std::endl;
      print_mat("%10.2e", m, l, lda, perm);
      std::cout << "D:" << std::endl;
      print_d<double>(m, d);
   }

   // Perform solve
   double *soln = new double[m];
   solve(m, q1, perm, l, lda, d, b, soln);
   if(debug) {
      printf("soln = ");
      for(int i=0; i<m; i++) printf(" %le", soln[i]);
      printf("\n");
   }

   // Check residual
   double bwderr = backward_error(m, a, lda, b, 1, soln, m);
   if(debug) printf("bwderr = %le\n", bwderr);
   EXPECT_LE(bwderr, 2e-13) << "(test " << test << " seed " << seed << ")" << std::endl;

   // Cleanup memory
   delete[] a; delete[] l;
   delete[] b;
   delete[] perm;
   delete[] work;
   delete[] d; delete[] soln;

   return failed ? -1 : 0;
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

int ldlt_tpp_torture_test(double u, double small, int ntest, int m, int n, bool debug=false) {
   for(int test=0; test<ntest; test++) {
      // Record seed we're using
      unsigned int seed = rand();
      //seed = 291189613;
      srand(seed);

      // Matrix has:
      //    70% chance of getting delays
      //    20% chance of being made singular
      bool delays = ( ((float) rand())/RAND_MAX < 0.7 );
      bool singular = ( ((float) rand())/RAND_MAX < 0.2 );

      // Output debug info
      if(debug) {
         std::cout << "##########################################" << std::endl;
         std::cout << "# Test " << test << " singular:" << (singular?"T":"F") << " delays:" << (delays?"T":"F") << std::endl;
         std::cout << "# Random seed " << seed << std::endl;
         std::cout << "##########################################" << std::endl;
      }

      int err = ldlt_tpp_test(u, small, delays, singular, m, n, debug, test, seed);
      if(err!=0) return err;
   }

   return 0; // Success
}

} /* anon namespace */

int run_ldlt_tpp_tests() {
   int nerr = 0;

   /* Simple tests */
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 1, 1) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 2, 2) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 3, 3) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 2, 1) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 3, 2) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 5, 3) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 8, 4) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, false, 33, 21) ));

   /* With delays */
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, false, 8, 4) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, false, 12, 3) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, false, 29, 7) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, false, 233, 122) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, false, 500, 500) ));

   /* Singular */
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, true, 8, 4) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, true, 12, 3) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, true, 29, 7) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, true, 233, 122) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, false, true, 500, 500) ));

   /* Singular, with delays */
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, true, 8, 4) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, true, 12, 3) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, true, 29, 7) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, true, 233, 122) ));
   TEST(( ldlt_tpp_test(0.01, 1e-20, true, true, 500, 500) ));

   /* Torture tests */
   TEST(( ldlt_tpp_torture_test(0.01, 1e-20, 1000, 100, 100) ));
   TEST(( ldlt_tpp_torture_test(0.01, 1e-20, 1000, 100, 50) ));

   return nerr;
}
