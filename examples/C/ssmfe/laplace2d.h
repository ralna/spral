#pragma once

/* Multiply x by Laplacian matrix for a square of size mx by my */
void apply_laplacian(
      int mx, int my, int nx, const double x[nx][my][mx], double Ax[nx][my][mx]
      ) {
   for(int k=0; k<nx; k++) {
      for(int i=0; i<mx; i++) {
         for(int j=0; j<my; j++) {
            double z = 4*x[k][j][i];
            if( i > 1  ) z -= x[k][j][i-1];
            if( j > 1  ) z -= x[k][j-1][i];
            if( i < mx ) z -= x[k][j][i+1];
            if( j < my ) z -= x[k][j+1][i];
            Ax[k][j][i] = z;
         }
      }
   }
}

/* Laplacian matrix for a rectangle of size nx x ny */
void set_laplacian_matrix(
      int nx, int ny, int lda, double a[nx*ny][lda]
      ) {
   for(int j=0; j<nx*ny; j++)
   for(int i=0; i<lda; i++)
      a[j][i] = 0.0;
   for(int ix=0; ix<nx; ix++) {
      for(int iy=0; iy<ny; iy++) {
        int i = ix + (iy - 1)*nx;
        a[i][i] = 4;
        if( ix >  1 ) a[i -  1][i] = -1;
        if( ix < nx ) a[i +  1][i] = -1;
        if( iy > 1  ) a[i - nx][i] = -1;
        if( iy < ny ) a[i + nx][i] = -1;
      }
   }
}

/* Apply one Gauss-Seidel step for preconditioning */
void apply_gauss_seidel_step(
      int mx, int my, int nx, const double x[nx][my][mx], double Tx[nx][my][mx]
      ) {
   for(int i=0; i<mx; i++)
      for(int j=0; j<my; j++)
         for(int k=0; k<nx; k++)
            Tx[k][j][i] = x[k][j][i]/4;
   for(int k=0; k<nx; k++) {
      /* forward update */
      for(int i=0; i<mx; i++) {
         for(int j=0; j<my; j++) {
            double z = 0.0;
            if( i > 1  ) z += Tx[k][j][i-1];
            if( j > 1  ) z += Tx[k][j-1][i];
            if( i < mx ) z += Tx[k][j][i+1];
            if( j < my ) z += Tx[k][j+1][i];
            Tx[k][j][i] = (Tx[k][j][i] + z/4)/4;
         }
      }
      /* backward update */
      for(int i=0; i<mx; i++) {
         for(int j=0; j<my; j++) {
            double z = 0.0;
            if( i > 1  ) z += Tx[k][j][i-1];
            if( j > 1  ) z += Tx[k][j-1][i];
            if( i < mx ) z += Tx[k][j][i+1];
            if( j < my ) z += Tx[k][j+1][i];
            Tx[k][j][i] = (Tx[k][j][i] + z/4)/4;
         }
      }
   }
}
