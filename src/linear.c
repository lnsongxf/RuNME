/* 
 * Solving Systems of Linear Equations
 * Runnan Yang 2013
 */

#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"

// Given an nxn triangular matrix A and nx1 vector B, solve Ax=B for nxn matrix x
//void BackSubIterate(double* A, double* B, int n, double** x){
//  
//}

void BackSubstitution(double* A, double* B, int n, double** x, int lt){
  
  // allocate x
  double* xtemp;
  zeros(n, 1, &xtemp);

  int sum = 0;

  if (lt){
    // If Lower Triangular

    for (int k = 0; k < n; k++){
      if (k == 0){
	xtemp[0] = B[0] / A[0*n + 0];
      }
      else{
	sum = 0;

	for (int j = 0; j < k - 1; j++){
	  sum += A[k*n + j] * xtemp[j];
	}

	xtemp[k] = (B[k] - sum) / A[k*n + k];
      }
      
    }

  }
  else{
    // If Upper Triangular
    for (int k = n-1; k >= 0; k--){
      if (k == n-1){
	xtemp[n-1] = B[n-1] / A[n-1];
      }
      else{
	sum = 0;

	for (int j = k+1; j < n; j++){
	  sum += A[k*n + j] * xtemp[j];
	}

	xtemp[k] = (B[k] - sum) / A[k*n + k];
      }
      
    }
  }
  copy(xtemp, x, n, 1);

}

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

  //printf("L\n");
  //print(Ltemp,2,2);

  // Get I - L
  double* diff;
  subtract(id, Ltemp, n, n, &diff);
  
  //printf("diff =id - lt\n");
  //print(diff, 2, 2);

  // A2 = (I - L) * A
  mult(diff, A, n, n, n, n, &A2temp);

  // Copy the temp matrices into output matrices
  //printf("now copying\n");
  copy(diff, L, n, n);
  //printf("...");
  copy(A2temp, A2, n, n);
  //printf("copying done\n");

  //printf("A\n");
  //print(A,2,2);

  //printf("A2\n");
  //print(A2temp,2,2);

  //printf("A2 = diff * A\n");
  //print(*A2,2,2);

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
      //printf("LU %d start\n", j);
      LUIterate(matrix, n, 0, &Ltemp, &A2);
      //printf("LU %d done\n", j);

      copy(Ltemp, &Lprod, n, n);
      //printf("Copy L %d done\n", j);
    }
    else{
      //printf("LU %d start\n", j);
      LUIterate(A, n, j, &Ltemp, &A2);
      //printf("LU %d done\n", j);

      mult(Lprod, Ltemp, n, n, n, n, &Lprodnew);
      //printf("Mult done\n");
      copy(Lprodnew, &Lprod, n, n);
      //printf("Copy L %d done\n", j);
    }

    // Update A
    copy(A2, &A, n, n);
    //printf("Copy A %d done\n", j);
    //printf("====================\n");

  }

  // Copy out
  copy(A, Uout, n, n);
  copy(Lprod, Lout, n, n);

}

void CholeskyDecomposition(double* A, int n, double** L){

}

void GaussJacobiIteration(double* A, double* b, double* x, int n, double** x2){

  double* x2temp;
  zeros(n, 1, &x2temp);
  double sum;

  // Gauss-Jacobi Iteration Formula
  // x^{k+1}_{j} = (1/A_{j,j}) (b_j - \sum_{i \neq j} (A_{j,i} * x_i))

  for (int i = 0; i < n; i++){
    sum = 0;

    for (int j = 0; j < n; j++){

      if (i != j){
	printf("[%f][%f]\n", A[i*n + j], x[j]);
	sum += A[i*n + j] * x[j];
      }
    }

    printf("%f\n", sum);

    x2temp[i] = (b[i] - sum) / A[i*n + i];
  }

  copy(x2temp, x2, n, n);
}

/*
*
* Test Cases
*
*/

void testLU(){
  double *mat;

  zeros(2,2,&mat);

  mat[0] = 1;
  mat[1] = 2;
  mat[1*2] = 3;
  mat[1*2+1] = 4;

  double* b;
  zeros(2,1,&b);

  b[0] = 1;
  b[1] = 2;

  double* L, *U, *z, *x;

  zeros(2,2,&L);
  zeros(2,2,&U);
  zeros(2,1,&z);
  zeros(2,1,&x);

  printf("A\n");
  print(mat,2,2);

  printf("B\n");
  print(b,2,1);

  LUDecomposition(mat, 2, &L, &U);

  printf("L\n");
  print(L,2,2);
  printf("U\n");
  print(U,2,2);

  // Solve Lz = b
  BackSubstitution(L, b, 2, &z, 1);

  // Solve Ux = z
  BackSubstitution(U, z, 2, &x, 0);

  printf("x\n");
  print(x,2,1);

  // Ax
  double* bstar;
  zeros(2,1,&bstar);

  mult(mat,x,2,2,2,1,&bstar);
  print(bstar,2,1);
}

void testGJ(){
  double *mat;

  zeros(2,2,&mat);

  mat[0] = 1;
  mat[1] = 2;
  mat[1*2] = 3;
  mat[1*2+1] = 4;

  double* b;
  zeros(2,1,&b);

  b[0] = 1;
  b[1] = 2;

  double* x0;
  zeros(2,1,&x0);

  double* x;
  zeros(2,1,&x);

  double* xtemp;
  zeros(2,1,&xtemp);

  printf("A\n");
  print(mat,2,2);

  printf("B\n");
  print(b,2,1);

  printf("x0\n");
  print(x0,2,1);

  int iters = 5;

  for (int i = 0; i < iters; i++){
    if (i == 0){
      GaussJacobiIteration(mat, b, x0, 2, &x);
    }
    else{
      GaussJacobiIteration(mat, b, x, 2, &xtemp);
      copy(xtemp,&x,2,1);
    }
    printf("Iteration %d:\n", i);
    print(x,2,1);
  }
}


int main(){
  testLU();
  //testGJ();

}
