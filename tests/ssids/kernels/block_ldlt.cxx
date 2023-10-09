/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * Licence: BSD licence, see LICENCE file for details
 *
 */
#include "block_ldlt.hxx"

#include <climits>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>

#include "framework.hxx"
#include "AlignedAllocator.hxx"
#include "ssids/cpu/kernels/wrappers.hxx"
#include "ssids/cpu/kernels/block_ldlt.hxx"

using namespace spral::ssids::cpu;
using namespace spral::test;

template<typename T>
void sym_mv(int n, const T *a, int lda, const T *x, T *b) {
   for(int i=0; i<n; i++) b[i] = 0.0;
   for(int c=0; c<n; c++) {
      b[c] += a[c*lda+c] * x[c];
      for(int r=c+1; r<n; r++) {
         b[c] += a[c*lda+r] * x[r];
         b[r] += a[c*lda+r] * x[c];
      }
   }
}

template <typename T>
void host_dslv(int n, const T *d, T *b) {
   for(int i=0; i<n; ) {
      if(d[2*i+1]==0.0) {
         // 1x1 pivot
         T d11 = d[2*i];
         b[i] *= d11;
         i++;
      } else {
         // 2x2 pivot
         T d11 = d[2*i];
         T d21 = d[2*i+1];
         T d22 = d[2*i+3];
         T b1 = b[i];
         T b2 = b[i+1];
         b[i]   = d11*b1 + d21*b2;
         b[i+1] = d21*b1 + d22*b2;
         i += 2;
      }
   }
}

template<typename T>
void solve(int n, const int *perm, const T *l, int ldl, const T *d, const T *b, T *x) {
   for(int i=0; i<n; i++) x[i] = b[perm[i]];
   // Fwd slv
   host_trsv<T>(FILL_MODE_LWR, OP_N, DIAG_UNIT, n, l, ldl, x, 1);
   // Diag slv
   host_dslv<T>(n, d, x);
   // Bwd slv
   host_trsv<T>(FILL_MODE_LWR, OP_T, DIAG_UNIT, n, l, ldl, x, 1);
   // Undo permutation
   T *temp = new T[n];
   for(int i=0; i<n; i++) temp[i] = x[i];
   for(int i=0; i<n; i++) x[perm[i]] = temp[i];
   // Free mem
   delete[] temp;
}

template <typename T>
T infnorm(int n, const T* a, int lda) {
   // || A ||_\infty = \max_{1\le i\le n} \sum_{j=1}^n | a_{ij} |
   T *rsum = new T[n];
   for(int i=0; i<n; i++) rsum[0] = 0.0;
   for(int j=0; j<n; j++) {
      for(int i=j; i<n; i++) {
         rsum[i] += fabs(a[j*lda+i]);
         if(i!=j) rsum[j] += fabs(a[j*lda+i]);
      }
   }
   T maxval=0.0;
   for(int i=0; i<n; i++)
      if(rsum[i] > maxval) maxval = rsum[i];
   delete[] rsum;
   return maxval;
}

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

template<typename T, int BLOCK_SIZE>
void find_maxloc_simple(const int from, const T *a, int lda, T &bestv, int &rloc, int &cloc, int &count) {
   bestv = -1.0;
   rloc = INT_MAX; cloc = INT_MAX;
   count = 0;
   for(int c=from; c<BLOCK_SIZE; c++) {
      for(int r=c; r<BLOCK_SIZE; r++) {
         T v = a[c*lda+r];
         if(fabs(v) > bestv) {
            bestv = fabs(v);
            rloc = r;
            cloc = c;
            count = 1;
         } else if (fabs(v) == bestv) {
            count += 1;
         }
      }
   }
   if(cloc < BLOCK_SIZE && rloc < BLOCK_SIZE)
      bestv = a[cloc*lda+rloc];
   else
      bestv = 0.0;
}

template<typename T, int BLOCK_SIZE, bool debug=false>
int test_maxloc(int from, bool zero=false) {
   int const lda = 2*BLOCK_SIZE;
#if defined(__AVX512F__)
   alignas(64)
#elif defined(__AVX__)
   alignas(32)
#else
   alignas(16)
#endif
     T a[BLOCK_SIZE*lda];

   /* Setup a random matrix. Entries in lwr triangle < 1.0, others = 100.0 */
   for(int j=0; j<from; j++)
   for(int i=0; i<BLOCK_SIZE; i++)
      a[j*lda+i] = 100.0;
   for(int j=from; j<BLOCK_SIZE; j++) {
      for(int i=0; i<j; i++)
         a[j*lda+i] = 100.0;
      for(int i=j; i<BLOCK_SIZE; i++)
         a[j*lda+i] = zero ? 0.0 : 2*((T) rand())/RAND_MAX -1;
   }

   /* Print matrix for debug */
   if(debug && BLOCK_SIZE < 16) {
      printf("A:\n");
      for(int i=0; i<BLOCK_SIZE; i++) {
         printf("%d:", i);
         for(int j=0; j<BLOCK_SIZE; j++)
            printf(" %e", a[j*lda+i]);
         printf("\n");
      }
   }

   /* Call both simple and avx versions and check they get the same answer */
   T mv1, mv2;
   int rloc1, cloc1, rloc2, cloc2, count;

   find_maxloc_simple<T, BLOCK_SIZE>(from, a, lda, mv1, rloc1, cloc1, count);
   block_ldlt_internal::find_maxloc<T, BLOCK_SIZE>(from, a, lda, mv2, rloc2, cloc2);

   // Compare them
   ASSERT_LE(fabs(mv1), 1.0); // If this is wrong, find_maxloc_simple is wrong
   if(count > 1) {
      // if there are multiple entries with the same absolute value as the
      // maximum, we cannot assume that both versions will find it at the same
      // location, since they essentially look through the matrix in a different
      // order. Only thing we can check for is the absolute values being equal.
      ASSERT_EQ(fabs(mv1), fabs(mv2));
   } else {
      ASSERT_EQ(mv1, mv2);
      ASSERT_EQ(rloc1, rloc2);
      ASSERT_EQ(cloc1, cloc2);
   }

   return 0; // Success
}

template<typename T, int BLOCK_SIZE, bool debug=false>
int test_maxloc_torture(int ntest) {
   for(int i=0; i<ntest; i++) {
      // Record seed we're using
      unsigned int seed = rand();
      srand(seed);

      int from = rand() % BLOCK_SIZE;

      // Do the test
      int err = test_maxloc<T, BLOCK_SIZE, debug>(from);
      if(err!=0) return err;
   }

   return 0; // Success
}

template <typename T,
          int BLOCK_SIZE,
          bool singular, // Make matrix singular
          bool debug // Switch on debugging output
          >
int ldlt_test_block(T u, T small) {
   bool failed = false;
   int n = BLOCK_SIZE;
   bool action = true; // Don't abort in presence of singularity

   // Generate test matrix
   int const lda = 2*BLOCK_SIZE;
   T *a = new T[n*lda];
   gen_sym_indef(n, a, lda);
   if(singular) {
      int col1 = (2 < n) ? 2 : 0;
      int col2 = ((n-1)-2 >= 0) ? (n-1)-2 : n-1;
      make_singular(n, col1, col2, a, lda);
   }

   // Generate a RHS based on x=1, b=Ax
   T b[BLOCK_SIZE];
   gen_rhs(BLOCK_SIZE, a, lda, b);

   // Print out matrices if requested
   if(debug) {
      std::cout << "A:" << std::endl;
      print_mat("%10.2e", n, a, lda);
   }

   // Factorize using main routine
   AlignedAllocator<T> Talloc;
   T *l = Talloc.allocate(n*lda);
   memcpy(l, a, n*lda*sizeof(T)); // Copy a to l
   int perm[BLOCK_SIZE];
   for(int i=0; i<n; i++) perm[i] = i;
   T d[2*BLOCK_SIZE];
   T ld[BLOCK_SIZE*BLOCK_SIZE];
   block_ldlt<T, BLOCK_SIZE>(0, perm, l, lda, d, ld, action, u, small);

   // Print out matrices if requested
   if(debug) {
      std::cout << "L:" << std::endl;
      print_mat("%10.2e", n, l, lda, perm);
      std::cout << "D:" << std::endl;
      print_d<T>(n, d);
   }

   // Perform solve
   T soln[BLOCK_SIZE];
   solve(n, perm, l, lda, d, b, soln);

   // Check residual
   T bwderr = backward_error(n, a, lda, b, 1, soln, BLOCK_SIZE);
   if(debug) printf("bwderr = %le\n", bwderr);
   EXPECT_LE(bwderr, 1e-14);

   // Cleanup memory
   delete[] a;
   Talloc.deallocate(l, n*lda);

   return (failed) ? -1 : 0;
}

template <typename T,
          int BLOCK_SIZE,
          int nblock,
          bool debug
          >
int ldlt_block_torture_test(T u, T small) {
   int n = BLOCK_SIZE;
   bool failed = false;
   bool action = true; // Don't abort in case of singularity

   // Generate nblock test matrices with entries as Unif(-10,10)
   T *a = new T[n*n*nblock];
   ASSERT_TRUE(a);
   for(int i=0; i<n*n*nblock; i++)
      a[i] = 20*((T) rand())/RAND_MAX - 10;
   // Each matrix has a 20% chance of being made singular
   int col1 = (2 < n) ? 2 : 0;
   int col2 = ((n-1)-2 >= 0) ? (n-1)-2 : n-1;
   for(int blk=0; blk<nblock; blk++)
      if( ((float) rand())/RAND_MAX < 0.2 ) {
         make_singular<T>(n, col1, col2, &a[n*n*blk], n);
      }

   // Generate RHSs based on x=1, b=Ax
   T x[BLOCK_SIZE];
   T *b = new T[n*nblock];
   for(int i=0; i<n; i++) x[i] = i+1;
   for(int blk=0; blk<nblock; blk++)
      sym_mv<T>(n, &a[n*n*blk], n, x, &b[n*blk]);

   // Factorize using main routine
   AlignedAllocator<T> Talloc;
   T *l = Talloc.allocate(n*n*nblock);
   T *ld = Talloc.allocate(n*n*nblock);
   ASSERT_TRUE(l);
   for(int i=0; i<n*n*nblock; i++) l[i] = a[i];
   int *perm = new int[n*nblock];
   for(int i=0; i<n*nblock; i++) perm[i] = i % n;
   T *d = new T[2*n*nblock];
   for(int blk=0; blk<nblock; blk++)
      block_ldlt<T, BLOCK_SIZE>(0, &perm[n*blk], &l[n*n*blk], BLOCK_SIZE, &d[2*n*blk], &ld[n*n*blk], action, u, small);

   // Check each system's residual
   T *soln = new T[n];
   for(int blk=0; blk<nblock; blk++) {
      solve(n, &perm[blk*n], &l[blk*n*n], n, &d[blk*2*n], &b[blk*n], soln);

      // Check residual
      T bwderr = backward_error(n, &a[blk*n*n], n, &b[blk*n], 1, soln, n);
      EXPECT_LE(bwderr, 1e-14);
   }

   // Free host memory
   delete[] a;
   Talloc.deallocate(l, n*n*nblock); Talloc.deallocate(ld, n*n*nblock);
   delete[] d;
   delete[] b; delete[] soln;
   delete[] perm;

   return (failed) ? -1 : 0;
}

int run_block_ldlt_tests() {
   int nerr = 0;

   /* Unit tests */
   TEST((test_maxloc<double, 128>(0, true)));
   /* Simple tests, block */
   TEST((ldlt_test_block<double, 2, false, false>(0.01, 1e-20)));
   TEST((ldlt_test_block<double, 4, false, false>(0.01, 1e-20)));
   TEST((ldlt_test_block<double, 8, false, false>(0.01, 1e-20)));
   TEST((ldlt_test_block<double, 16, false, false>(0.01, 1e-20)));

   /* Singular matrix tests, block */
   TEST((ldlt_test_block<double, 2, true, false>(0.01, 1e-20)));
   TEST((ldlt_test_block<double, 4, true, false>(0.01, 1e-20)));
   TEST((ldlt_test_block<double, 8, true, false>(0.01, 1e-20)));
   TEST((ldlt_test_block<double, 16, true, false>(0.01, 1e-20)));


   /************************************************************
    * Torture Tests
    ************************************************************/
   TEST((test_maxloc_torture<double, 128>(10000)));
   TEST((ldlt_block_torture_test<double, 16, 500, false>(0.01, 1e-20)));

   return nerr;
}
