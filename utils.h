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
Matrix* get_matrix_from_2D_array (double** points, int num_points, int num_coordinates);
void print_matrix (Matrix* m);
void print_full_matrix (Matrix* m);
Matrix* get_zeros_matrix_size_n (int n);
Matrix* get_n_k_zero_matrix (int n, int k);
Matrix* get_identity_matrix_size_n (int n);
void get_cell_with_largest_value (Matrix* m, Cell* cell_pointer);
int is_diagonal_matrix (Matrix* m);
double sign (double num);
Matrix* multiply_matrices_and_free_memory (Matrix* m1, Matrix* m2);
void free_matrix_memory (Matrix* m);
Matrix* get_copy_of_matrix(Matrix* original);
Matrix* normalize_matrix (Matrix* original);
int get_num_points(FILE *fptr);
int get_num_coordinates(FILE *fptr);
void get_points_from_file (int num_points, int num_coordinates, FILE *fptr, double** points);
void single_line_to_point (double* point, char* single_line);

#endif //SPK_SPKMEANSMODULE_H