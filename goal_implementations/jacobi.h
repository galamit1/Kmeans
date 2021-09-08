#ifndef SPK_JACOBI_H
#define SPK_JACOBI_H

#include "../utils.c"
#include <assert.h>
#include <math.h>

/**********************/
/**FUNCTION DECLARATION**/
/**********************/
Matrix* run_jacobi (Matrix* lnorm);
void get_column_by_index(Matrix* m1, Matrix* m2, int m1_index, int m2_index, int column_size);
int get_eigen_gap_k (double* eigen_valus_array, int* indexes_array, int array_length);
void swap_doubles(double *xp, double *yp);
void swap_ints(int *xp, int *yp);
void bubble_sort_preserve_indexs(double* values, int* indexes, int array_length);
int* get_initial_indexes_array (int n);
Matrix* get_eigen_vectors_and_values (Matrix* original_matrix, double* eigen_valus_array);
Matrix* preform_pivot_step_and_free_memory (Matrix* A, int i, int j, double c, double s);
double get_off_diagonal_sum_of_matrix (Matrix* m);
Matrix* get_rotation_matrix_for_m (Matrix* m, Cell* largest_non_diagonal_cell, double c, double s);

#endif //SPK_JACOBI_H