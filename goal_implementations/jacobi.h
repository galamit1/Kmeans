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
void get_column_by_index(Matrix* m1, Matrix* m2, int m1_index, int m2_index, int column_size);
int getEigenGapK (double* eigen_valus_array, int* indexes_array, int array_length);
void swap_doubles(double *xp, double *yp);
void swap_ints(int *xp, int *yp);
void bubble_sort_preserve_indexs(double* values, int* indexes, int array_length);
int* get_initial_indexes_array (int n);
Matrix* getEigenVectorsAndValues (Matrix* originalMatrix, double* eigenValusArray);
Matrix* preformPivotStepAndFreeMemory (Matrix* A, int i, int j, double c, double s);
double getOffDiagonalSumOfMatrix (Matrix* m);
Matrix* getRotationMatrixForM (Matrix* m, Cell* largestNonDiagonalCell, double c, double s);