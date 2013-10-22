/* 
 * Solving Systems of Linear Equations
 * Runnan Yang 2013
 */

#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"

void LUIterate(double* A, int n, int curr, double** L, double** A2){

  // A[k+1] = (I - L) A[k]

  // Create identity matrix
  double* id;
  eye(n, &id);
  
  double* Ltemp;
  double* A2temp;

  zeros(n, n, &Ltemp);
  zeros(n, n, &A2temp);
  //zeros(n, n, L);
  
  // Set up L matrix
  for (int i = 0; i < n; i++){
    Ltemp[i*n + curr] = A[i*n + curr] / A[curr*n + curr];
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

  L = &Ltemp;
  A2 = &A2temp;

  printf("A\n");
  print(A,2,2);

  printf("A2\n");
  print(A2temp,2,2);

  printf("A2 = diff * A\n");
  print(*A2,2,2);

}

void LUDecomposition(double** matrix, int rows, int cols){
  
}

int main(){

  double *mat;

  zeros(2,2,&mat);

  mat[0] = 1;
  mat[1] = 2;
  mat[1*2] = 3;
  mat[1*2+1] = 4;

  //LUDecomposition(mat, 2, 2);

  double** L;
  double* A2;

  print(mat,2,2);

  LUIterate(mat, 2, 0, L, &A2);
  
  printf(">>>\t%f\n",);
}
