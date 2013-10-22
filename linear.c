/* 
 * Solving Systems of Linear Equations
 * Runnan Yang 2013
 */

#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"

void LUIterate(double* A, int n, int k, double** L, double** A2){
  // A[k+1] = (I - L) A[k]

  // Create identity matrix
  double* id;
  eye(n, &id);
  
  double* Ltemp;
  double* A2temp;

  zeros(n, n, &Ltemp);
  zeros(n, n, &A2temp);
  
  // Set up L matrix
  for (int i = k + 1; i < n; i++){
    Ltemp[i*n + k] = A[i*n + k] / A[k*n + k];
  }

  printf("L\n");
  print(Ltemp,2,2);

  // Get I - L
  double* diff;
  subtract(id, Ltemp, n, n, &diff);
  
  printf("diff =id - lt\n");
  print(diff, 2, 2);

  // A2 = (I - L) * A
  mult(diff, A, n, n, n, n, &A2temp);

  // Copy the temp matrices into output matrices
  printf("now copying\n");
  copy(diff, L, n, n);
  printf("...");
  copy(A2temp, A2, n, n);
  printf("copying done\n");

  printf("A\n");
  print(A,2,2);

  printf("A2\n");
  print(A2temp,2,2);

  printf("A2 = diff * A\n");
  print(*A2,2,2);

}

void LUDecomposition(double* matrix, int n, double** Lout, double** Uout){
  // for each column

  double* Lprod;
  double* Ltemp;

  double* Lprodnew;

  double* A;
  double* A2;

  zeros(n, n, &Ltemp);
  zeros(n, n, &Lprod);
  zeros(n, n, &A);
  zeros(n, n, &A2);
  zeros(n, n, &Lprodnew);
  
  for (int j = 0; j < n - 1; j++){

    if (j == 0){
      printf("LU %d start\n", j);
      LUIterate(matrix, n, 0, &Ltemp, &A2);
      printf("LU %d done\n", j);

      copy(Ltemp, &Lprod, n, n);
      printf("Copy L %d done\n", j);
    }
    else{
      printf("LU %d start\n", j);
      LUIterate(A, n, j, &Ltemp, &A2);
      printf("LU %d done\n", j);

      mult(Lprod, Ltemp, n, n, n, n, &Lprodnew);
      printf("Mult done\n");
      copy(Lprodnew, &Lprod, n, n);
      printf("Copy L %d done\n", j);
    }

    // Update A
    copy(A2, &A, n, n);
    printf("Copy A %d done\n", j);
    printf("====================\n");

  }

  // Copy out
  copy(A, Uout, n, n);
  copy(Lprod, Lout, n, n);

}

int main(){

  double *mat;

  zeros(2,2,&mat);

  mat[0] = 1;
  mat[1] = 2;
  mat[1*2] = 3;
  mat[1*2+1] = 4;

  //LUDecomposition(mat, 2, 2);

  double* L;
  double* U;

  zeros(2,2,&L);
  zeros(2,2,&U);

  print(mat,2,2);

  //LUIterate(mat, 2, 0, &L, &A2);

  LUDecomposition(mat, 2, &L, &U);

  print(L,2,2);
  print(U,2,2);
}
