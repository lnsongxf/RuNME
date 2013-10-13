int add(double** matrixA, double** matrixB, int rows, int cols, double*** sum);

int subtract(double** matrixA, double** matrixB, int rows, int cols, double*** diff);

int multiply(double** matrixA, double** matrixB, int rowsA, int colsA, int rowsB, int colsB, double*** prod);

int ones(int n, int m, double*** one);
int zeros(int n, int m, double*** zero);
int eye(int n, double*** id);

int print(int n, int m, double** matrix);
