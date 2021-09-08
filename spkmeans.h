#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define LINE_SIZE 1001 /*Given assumption each line is no more than 1000 characters + \n character*/
#define MAXIMUM_DOUBLE 1.7E+308
#define MAX_ITER 300
#define EPSILON 1.0e-15
#define NUM_ITERATIONS 100

/**********************/
/**STRUCT DEFINITIONS**/
/**********************/
typedef struct cluster {
    int cluster_size;
    int cluster_index;
    double* prev_centroids;
    double* curr_centroids;
    double* sum;
} Cluster;

typedef struct Matrix {
    double** cells;
    int rows;
    int cols;
} Matrix ;

typedef struct CellStruct {
    int row;
    int col;
    double value;
} Cell;

/**********************/
/**FUNCTION DECLARATION**/
/**********************/

int main(int argc, char **argv);
void run_functions_according_to_goal(char * goal, Matrix * points_matrix, int k);
/***UTILS***/
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
void multiply_matrices_to_existing_pointer (Matrix* m1, Matrix* m2, Matrix* product);
void free_matrix_memory (Matrix* m);
Matrix* get_copy_of_matrix(Matrix* original);
Matrix* normalize_matrix (Matrix* original);
int get_num_points(FILE *fptr);
int get_num_coordinates(FILE *fptr);
void get_points_from_file (FILE *fptr, double** points);
void single_line_to_point (double* point, char* single_line);
/*** WAM ***/
Matrix * run_wam(Matrix * points);
double calculate_weight(double * point1, double * point2, int num_coordinates);
/*** DDG ***/
Matrix * run_ddg(Matrix * wam);
void convert_ddg_with_the_pow_of_minus_half(Matrix * ddg_matrix);
/*** LNORM ***/
Matrix * run_lnorm(Matrix * wam, Matrix * ddg);
/*** JACOBI ***/
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
/*** SPK ***/
void run_spk(Matrix * points, int k);
Cluster** init_k_clusters (Matrix * points, int k);
double** init_points (int num_points, int num_coordinates);
double* convert_line_to_point (double* point, char* line);
Cluster* make_cluster (const double* point, const int num_coordinates, const int index);
void add_point_to_cluster (Cluster* cluster, const double* point, int num_coordinates);
void find_minimum_distance (Cluster** clusters, double* point, int k, int num_coordinates);
double update_cluster_centroids_and_sum (Cluster* cluster, int num_coordinates);
int update_all_clusters (Cluster** clusters, int k, int num_coordinates);
void free_clusters_memory (Cluster** clusters, int k);
void free_points_memory (double** points, int num_points);
double get_distance (Cluster* cluster, const double* point, int num_coordinates);