#include<stdio.h>
#include<stdlib.h>

int add(double* A, double* B, int rows, int cols, double** sum){
  *sum = calloc(rows * cols, sizeof(double));

  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      (*sum)[i*cols + j] = A[i*cols + j] + B[i*cols+j];
    }
  }

  return 0;
}

int subtract(double* A, double* B, int rows, int cols, double** diff){
  *diff = calloc(rows * cols, sizeof(double));

  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      (*diff)[i*cols + j] = A[i*cols + j] - B[i*cols+j];
    }
  }

  return 0;
}

int mult(double* A, double* B, int rowsA, int colsA, int rowsB, int colsB, double** prod){

  *prod = calloc(rowsA * colsB, sizeof(double));
  //printf("%d x %d\n", rowsA, colsB);
  // Check that parameters are sane!
  if (colsB < 1)
    return 0;
  // for each row of A
  for (int i = 0; i < rowsA; i++){

    // for each col of B
    for (int j = 0; j < colsB; j++){
      //printf("[%d, %d] = ", i, j);

      double sum = 0;

      // Dot product of ith row of A, jth col of B
      // for the kth element
      for (int k = 0; k < colsA; k++){
	// sum = A[i,k] + B[k,j]
	sum += A[i*colsA + k] * B[k*colsB + j];
      }

      // AB[i,j] = sum
      (*prod)[i*colsB + j] = sum;
      //printf("%f\n", sum);
    }
  }

  return 0;
}

int ones(int n, int m, double** one){
  *one = calloc(n*m, sizeof(double));

  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      (*one)[i*m + n] = 1;
    }
  }

  return 0;
}

int zeros(int n, int m, double** zero){
  *zero = calloc(n*m, sizeof(double));
  /*
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      (*zero)[i*m + n] = 0;
    }
  }
  */
  return 0;
}

int eye(int n, double** id){
  *id = calloc(n*n, sizeof(double));

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i == j){
	(*id)[i*n + j] = 1;
      }
      else{
	(*id)[i*n + j] = 0;
      }
    }
  }

  return 0;
}

int print(double* A, int r, int c){
  for (int i = 0; i < r; i++){
    for (int j = 0; j < c; j++){
      printf("%f\t",A[i*c + j]);
    }
    printf("\n");
  }
  printf("\n");

  return 0;
}

/*
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
*/

/*
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
*/
/*
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
*/


/*
int main(){
  int n = 2;
  int m = 2;

  double* x;
  double* y;

  ones(2,2,&x);
  eye(2,&y);

  x[0] = 5;

  double* z;

  add(x,y,2,2,&z);
  mult(x,y,2,2,2,2,&z);

  print(x,2,2);
  print(y,2,2);

  print(z,2,2);
}
*/
