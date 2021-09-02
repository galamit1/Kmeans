#ifndef SPK_SPKMEANSMODULE_H
#define SPK_SPKMEANSMODULE_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

#endif //SPK_SPKMEANSMODULE_H