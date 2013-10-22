
int copy(double* matrixA, double** matrixB, int rows, int cols);

int add(double* matrixA, double* matrixB, int rows, int cols, double** sum);
int subtract(double* matrixA, double* matrixB, int rows, int cols, double** diff);
int mult(double* matrixA, double* matrixB, int rowsA, int colsA, int rowsB, int colsB, double** prod);

int ones(int n, int m, double** one);
int zeros(int n, int m, double** zero);
int eye(int n, double** id);

int print(double* matrix, int n, int m);
