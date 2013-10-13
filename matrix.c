#include<stdio.h>
#include<stdlib.h>

int add(double** matrixA, double** matrixB, int rows, int cols, double*** sum){
  
  // Create sum matrix
  
  *sum = calloc(rows, sizeof(double*));
  for (int i = 0; i < rows; i++){
    (*sum)[i] = (double*) calloc(cols, sizeof(double));
  }
  

  // Assign entries of matrix
  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      (*sum)[i][j] = matrixA[i][j] + matrixB[i][j];
    }
  }

  return 0;
}

int subtract(double** matrixA, double** matrixB, int rows, int cols, double*** diff){
  
  // Create diff matrix
  
  *diff = calloc(rows, sizeof(double*));
  for (int i = 0; i < rows; i++){
    (*diff)[i] = (double*) calloc(cols, sizeof(double));
  }
  

  // Assign entries of matrix
  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      (*diff)[i][j] = matrixA[i][j] - matrixB[i][j];
    }
  }

  return 0;
}

int multiply(double** matrixA, double** matrixB, int rowsA, int colsA, int rowsB, int colsB, double*** prod){
  
  // Check if valid multiplication
  if (colsA != rowsB)
    return 1;

  // Create product matrix
  *prod = calloc(rowsA, sizeof(double*));
  for (int i = 0; i < rowsA; i++){
    (*prod)[i] = (double*) calloc(colsB, sizeof(double));
  }

  // Assign entries of matrix
  for (int i = 0; i < rowsA; i++){
    for (int j = 0; j < colsB; j++){
      
      // row i dot col j
      int sum = 0;
      for (int k = 0; k < colsA; k++){
	sum += matrixA[i][k] * matrixB[k][j];
      }

      (*prod)[i][j] = sum;

    }
  }
  
  return 0;
}

// Ones Matrix
int ones(int n, int m, double*** one){

  *one = (double**) malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++){
    (*one)[i] = (double*) malloc(m * sizeof(double));
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      (*one)[i][j] = 1;
    }
  }
    
  return 0;
}

// Zero Matrix
int zeros(int n, int m, double*** zero){

  *zero = (double**) malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++){
    (*zero)[i] = (double*) malloc(m * sizeof(double));
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      (*zero)[i][j] = 0;
    }
  }
    
  return 0;
}

// Identity Matrix
int eye(int n, double*** id){

  *id = (double**) malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++){
    (*id)[i] = (double*) malloc(n * sizeof(double));
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i == j){
	(*id)[i][j] = 1;
      }
      else{
	(*id)[i][j] = 0;
      }
    }
  }

  return 0;
}

// Print 2d Matrix
int print(int n, int m, double** matrix){
  printf("\n");
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      printf("%f\t", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  return 0;
}
/*
int main(){
  int n = 2;
  int m = 2;

  double** x;
  double** y;

  int i = ones(2,2,&x);
  int j = eye(2,&y);

  x[0][0] = 5;

  double** z;

  //int k = add(x,y,2,2,&z);
  int k = multiply(x,y,2,2,2,2,&z);

  print(2,2,x);
  print(2,2,y);

  print(2,2,z);
}
*/
