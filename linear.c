/* 
 * Solving Systems of Linear Equations
 * Runnan Yang 2013
 */

#include <stdlib.h>
#include <stdio.h>

void LUDecomposition(double** matrix, int rows, int cols){
  // Just print for now
  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j ++){
      printf("%d", matrix[i][j]);
    }
  }

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

  LUDecomposition(mat, 2, 2);


  printf("hi");
}
