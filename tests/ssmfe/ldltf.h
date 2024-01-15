/* C Wrappers for LAPACK routines */
void dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
int cwrap_dsytrf(char uplo, int n, double *a, int lda, int *ipiv, double *work, int lwork) {
   int info;
   dsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
   return info;
}
void dsytrs_(char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
int cwrap_dsytrs(char uplo, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {
   int info;
   dsytrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
   return info;
}

/* Counts negative eigenvalues of factor D */
int num_neg_D(int n, int ld, double (*LDLT)[n][ld], int *ipiv) {
   int nneg = 0;
   for(int i=0; i<n; ) {
      double s = (*LDLT)[i][i];
      if( ipiv[i] < 0 ) {
         double t = (*LDLT)[i][i+1];
         double r = (*LDLT)[i+1][i+1];
         if( s*r - t*t < 0.0 ) {
            nneg++;
         } else if ( s*r - t*t > 0.0 && s + r < 0.0 ){
            nneg+=2;
         }
         i += 2;
      } else {
         if( s < 0.0 ) nneg++;
         i++;
      }
   }
   return nneg;
}

int num_neg_Dz(int n, int ld, double complex (*LDLT)[n][ld], int *ipiv) {
   int nneg = 0;
   for(int i=0; i<n; ) {
      double s = (*LDLT)[i][i];
      if( ipiv[i] < 0 ) {
         double t = fabs((*LDLT)[i][i+1]);
         double r = (*LDLT)[i+1][i+1];
         if( s*r - t*t < 0.0 ) {
            nneg++;
         } else if ( s*r - t*t > 0.0 && s + r < 0.0 ){
            nneg+=2;
         }
         i += 2;
      } else {
         if( s < 0.0 ) nneg++;
         i++;
      }
   }
   return nneg;
}

