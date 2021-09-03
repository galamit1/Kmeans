//
// Created by galam on 26/07/2021.
//

#include "../utils.c"
#include <assert.h>
#include <math.h>

/**********************/
/**FUNCTION DECLARATION**/
/**********************/
Matrix* run_jacobi (Matrix* lnorm);
int getEigenGapK (double* eigen_valus_array, int array_length);
int standart_sort (double a, double b);
Matrix* getEigenVectorsAndValues (Matrix* originalMatrix, double* eigenValusArray);
Matrix* preformPivotStepAndFreeMemory (Matrix* A, int i, int j, double c, double s);
double getOffDiagonalSumOfMatrix (Matrix* m);
Matrix* getRotationMatrixForM (Matrix* m, Cell* largestNonDiagonalCell, double c, double s);