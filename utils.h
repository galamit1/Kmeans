#ifndef SPK_SPKMEANSMODULE_H
#define SPK_SPKMEANSMODULE_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXIMUM_DOUBLE 1.7E+308

/**********************/
/**STRUCT DEFINITIONS**/
/**********************/


 struct MatrixStruct {
    double** cells;
    int rows;
    int cols;
}typedef Matrix ;

typedef struct CellStruct {
    int row;
    int col;
    double value;
} Cell;

/**********************/
/**FUNCTION DECLARATION**/
/**********************/

Matrix* get_points_matrix (char *path);
Matrix* getMatrixFrom2DArray (double** points, int numPoints, int numCoordinates);
void printMatrix (Matrix* m);
void printFullMatrix (Matrix* m);
Matrix* getZerosMatrixSizeN (int n);
Matrix* getIdentitiyMatrixSizeN (int n);
void getCellWithLargestValue (Matrix* m, Cell* cell_pointer);
int isDiagonalMatrix (Matrix* m);
double sign (double num);
Matrix* multiplyMatricesAndFreeMemory (Matrix* m1, Matrix* m2);
void freeMatrixMemory (Matrix* m);
Matrix* get_copy_of_matrix(Matrix* original);
int getNumPoints(FILE *fptr);
int getNumCoordinates(FILE *fptr);
void getPointsFromFile (int numPoints, int numCoordinates, FILE *fptr, double** points);
void singleLineToPoint (double* point, char* singleLine);

#endif //SPK_SPKMEANSMODULE_H