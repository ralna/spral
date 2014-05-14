#include "stdio.h"
#include "stdlib.h"
#include "cublas_v2.h"

#include "src/cuda/cuda_check.h"

float tdiff(struct timespec tp1, struct timespec tp2) {
   return (tp2.tv_sec - tp1.tv_sec) + 1e-9*(tp2.tv_nsec - tp1.tv_nsec);
}

int bench_nostream(int maxdimn, double *gpu_work) {
   char buffer[10000];
   cublasHandle_t cublas;
   cublasCreate(&cublas);

   struct timespec start, end;
   cudaDeviceSynchronize();
   clock_gettime(CLOCK_REALTIME, &start);

   long flops = 0;
   while(!feof(stdin)) {
      fgets(buffer, sizeof buffer, stdin);
      if(buffer[0] == '#') continue; // Skip comments
      // Read data
      int node, lev, blkm, blkn;
      sscanf(buffer, "%d %d %d %d", &node, &lev, &blkm, &blkn);
      long m = blkm-blkn;
      if(2*m*blkn+m*m > maxdimn) {
         printf("Can't cope with node %d being %d x %d\n", node, blkm, blkn);
         return 1;
      }
      // Launch BLAS
      double alpha = -1.0;
      double beta = 1.0;
      double *a = gpu_work;
      double *b = a + m*blkn;
      double *c = b + m*blkn;
      if(blkm-blkn != 0 && blkn != 0)
         CublasSafeCall(
            cublasDgemm(cublas, CUBLAS_OP_N, CUBLAS_OP_T, blkm-blkn, blkm-blkn,
               blkn, &alpha, a, blkm-blkn, b, blkm-blkn, &beta, c, blkm-blkn)
            );
      flops += 2*m*m*blkn;
   }
   cudaDeviceSynchronize();
   clock_gettime(CLOCK_REALTIME, &end);
   printf("nostream: Did %f flops\n", (float) flops);
   printf("nostream: Took %f seconds\n", tdiff(start, end));
   printf("nostream: Rate %f Gflop/s\n", flops/(1e9*tdiff(start,end)));

   cublasDestroy(cublas);
   return 0;
}

int bench_stream(int maxdimn, double *gpu_work) {
   char buffer[10000];
   cudaStream_t stream[16];
   for(int i=0; i<16; i++) CudaSafeCall( cudaStreamCreate(&stream[i]) );
   cublasHandle_t cublas[16];
   for(int i=0; i<16; i++) CublasSafeCall( cublasCreate(&cublas[i]) );
   for(int i=0; i<16; i++) CublasSafeCall( cublasSetStream(cublas[i], stream[i]) );

   struct timespec start, end;
   cudaDeviceSynchronize();
   clock_gettime(CLOCK_REALTIME, &start);

   long flops = 0;
   int lastlev = 1;
   int s=0;
   while(!feof(stdin)) {
      fgets(buffer, sizeof buffer, stdin);
      if(buffer[0] == '#') continue; // Skip comments
      // Read data
      int node, lev, blkm, blkn;
      sscanf(buffer, "%d %d %d %d", &node, &lev, &blkm, &blkn);
      long m = blkm-blkn;
      if(2*m*blkn+m*m > maxdimn) {
         printf("Can't cope with node %d being %d x %d\n", node, blkm, blkn);
         return 1;
      }
      // Check for level change
      if(lev != lastlev) {
         cudaDeviceSynchronize();
         lastlev = lev;
      }
      // Launch BLAS
      double alpha = -1.0;
      double beta = 1.0;
      double *a = gpu_work;
      double *b = a + m*blkn;
      double *c = b + m*blkn;
      //printf("running on stream %d: %d %d %d\n", s, blkm-blkn, blkm-blkn, blkn);
      if(blkm-blkn != 0 && blkn != 0)
         CublasSafeCall(
            cublasDgemm(cublas[s], CUBLAS_OP_N, CUBLAS_OP_T, blkm-blkn,
               blkm-blkn, blkn, &alpha, a, blkm-blkn, b, blkm-blkn, &beta, c,
               blkm-blkn)
            );
      flops += 2*m*m*blkn;
      if(++s >= 16) s=0;
   }
   CudaSafeCall( cudaDeviceSynchronize() );
   clock_gettime(CLOCK_REALTIME, &end);
   printf("stream: Did %f flops\n", (float) flops);
   printf("stream: Took %f seconds\n", tdiff(start, end));
   printf("stream: Rate %f Gflop/s\n", flops/(1e9*tdiff(start,end)));

   for(int i=0; i<16; i++) {
      CublasSafeCall( cublasDestroy(cublas[i]) );
      CudaSafeCall( cudaStreamDestroy(stream[i]) );
   }
   return 0;
}

int main(void) {

   const int maxdimn = 3*5000*5000;

   double *work = (double *) malloc(maxdimn*sizeof(double));
   double *gpu_work;
   int err = cudaMalloc(&gpu_work, maxdimn*sizeof(double));
   if(err!=0) {
      printf("Failed to alloc\n");
      return 1;
   }

   srand((unsigned)time(NULL));
   for(int i=0; i<maxdimn; i++) work[i] = ((double) rand()) / RAND_MAX;

   cudaMemcpy(gpu_work, work, maxdimn*sizeof(double), cudaMemcpyHostToDevice);

   //if(!bench_nostream(maxdimn, gpu_work)) { return 1; }
   if(!bench_stream(maxdimn, gpu_work)) { return 1; }

   cudaFree(gpu_work);
   free(work);

   return 0;
}
