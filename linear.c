/* 
 * Solving Systems of Linear Equations
 * Runnan Yang 2013
 */

#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"

void LUIterate(double** A, int n, int curr, double*** L, double*** A2){

  // Create identity matrix
  double** id;
  eye(n, &id);
  
  zeros(n, n, L);
  zeros(n, n, A2);

  for (int i = 0; i < n; i++){
    (*L)[i][n] = A[i][n] / A[n][n];
  }

  double** diff;
  subtract(id, L, n, n, &diff);
  
  double** prod;
  multiply(diff, A, n,n,n,n, &prod);

}

void LUDecomposition(double** matrix, int rows, int cols){
  
}

int main(){

  double **mat;
  mat = (double **) malloc(2 * sizeof(double *));

  for (int i = 0; i < 2; i++){
    mat[i] = (double *) malloc(2 * sizeof(double));
  }
  
  mat[0][0] = 1;
  mat[0][1] = 2;
  mat[1][0] = 3;
  mat[1][1] = 4;

  //LUDecomposition(mat, 2, 2);

  LUIterate(mat, 2, 2);
  
  printf("hi");
}
