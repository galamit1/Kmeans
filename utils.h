#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**********************/
/**STRUCT DEFINITIONS**/
/**********************/

typedef struct Matrix {
    double** cells;
    int rows;
    int cols;
} Matrix ;

typedef struct Cell {
    int row;
    int col;
    double value;
} Cell;

/**********************/
/**FUNCTION DECLARATION**/
/**********************/

double** get_points (char *path);
Matrix* getMatrixFrom2DArray (double** points, int numPoints, int numCoordinates);
void printMatrix (Matrix* m);
Matrix* getZerosMatrixSizeN (int n);
Matrix* getIdentitiyMatrixSizeN (int n);
Cell* getCellWithLargestValue (Matrix* m);
int isDiagonalMatrix (Matrix* m);
int sign (double num);
Matrix* multiplyMatricesAndFreeMemory (Matrix* m1, Matrix* m2);
void freeMatrixMemory (Matrix* m);
int getNumPoints(FILE *fptr);
int getNumCoordinates(FILE *fptr);
void getPointsFromFile (int numPoints, int numCoordinates, FILE *fptr, double** points);
void singleLineToPoint (double* point, char* singleLine);
double** get_points (char *path);