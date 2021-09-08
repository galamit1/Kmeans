//
// Created by galam on 26/07/2021.
//

#include "./jacobi.h"

#define EPSILON 1.0e-15
#define NUM_ITERATIONS 20

/**********************/
/**TESTING ONLY - DON'T FORGET TO REMOVE**/
/**********************/

 int main () {
     Matrix* lnorm;
     Matrix* jacobi;
     Matrix* T;

     lnorm = get_points_matrix("./lnorm-example-after-1-iteration.txt");

    jacobi = run_jacobi(lnorm);
    printf("*********************** \n");
    printf("Final Jacobi U matrix is: \n");
    print_matrix(jacobi);
    T = normalize_matrix(jacobi);
    printf("*********************** \n");
    printf("Final T matrix is: \n");
    print_matrix(T);
    printf("Finished printing final matrix\n");

    free_matrix_memory(lnorm);
    printf("Finished freeing lnorm");
    free_matrix_memory(jacobi);
    printf("Finished freeing Jacobi");
    free_matrix_memory(T);

    printf("Bye Bye\n");
    return 0;
 }

/**********************/
/**TESTING ONLY - DON'T FORGET TO REMOVE**/
/**********************/

/*Recieves lnorm symmetric matrix, runs jacobi algorithm and returns n*k matrix U  containing the first k eigenvectors u1, . . . , uk of Lnorm columns. K is inferred by matrix size*/
Matrix* run_jacobi (Matrix* lnorm) {
    Matrix* eigen_vectors_matrix;
    Matrix* final_U_matrix;
    double* eigen_valus_array;
    int* indexes_array;
    int n;
    int k;
    int i,j;

    n = lnorm -> rows;
    /*allocate memory for eigenvalues array*/
    eigen_valus_array = (double*)malloc(n * sizeof(double));
    assert (eigen_valus_array != NULL); //TODO: add printf off error message
    /*get eigenvectors and update eigenvalues array accordingly*/
    eigen_vectors_matrix = get_eigen_vectors_and_values(lnorm, eigen_valus_array);
    /*get initial indexes array i.e. [0,1,...,n-1] */
    indexes_array = get_initial_indexes_array(n);
    /*get eigengap k and preserve original indexes during sort*/
    k = get_eigen_gap_k(eigen_valus_array, indexes_array, n);
    /*allocate memory for U - a n*k Matrix*/
    final_U_matrix = get_n_k_zero_matrix(n, k);

    /*point to k eigenvectors that correspond to k smallest eigenvalues*/
    for (i=0; i < k; i++) {
        get_column_by_index(final_U_matrix, eigen_vectors_matrix, i, indexes_array[i], n);
    }

    return final_U_matrix;
}

/*Receive 2 Matrixes, column size and column indexes for both matrixes, copy column from relevant index in m2 to relevant index in m1*/
void get_column_by_index(Matrix* m1, Matrix* m2, int m1_index, int m2_index, int column_size) {
    int i;

    for (i=0; i < column_size; i++) {
        m1 -> cells[i][m1_index] = m2 -> cells[i][m2_index];
    }
}

/*Recieve eigenvalues array & initial indexes array, calculate k by eigengap heuristic and update indexes array along with sorting values array*/
int get_eigen_gap_k (double* eigen_valus_array, int* indexes_array, int array_length) {
    double max_gap, current_gap;
    int max_index, current_index;

    /*Sort the array and preserve indexes*/
    bubble_sort_preserve_indexs(eigen_valus_array, indexes_array, array_length);

    /*Set first index to be initial max gap & index*/
    max_index = 0;
    max_gap = fabs(eigen_valus_array[0] - eigen_valus_array[1]);

    /*Iterate rest of eigenValues*/
    for (current_index = 0; current_index < (array_length - 1); current_index++) {
        current_gap = fabs(eigen_valus_array[current_index] - eigen_valus_array[current_index + 1]);
        /*In case of equality in the argmax of some eigengaps, use the lowest index - therefore we use > and not >=*/
        if (current_gap > max_gap) {
            max_gap = current_gap;
            max_index = current_index;
        }
    }

    /*Indexes start with 0 but we need number of eigenvalues before gap - therefore we add 1 to max_index*/
    return (max_index + 1);
}


/*Bubble sort code from: https://www.geeksforgeeks.org/bubble-sort/*/
/*Recieve 2 double pointers and swap them*/
void swap_doubles(double *xp, double *yp) {
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void swap_ints(int *xp, int *yp) {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

/*Recieves an array of doubles & array of indexes, bubble-sorts both arrays accordingly in-place*/
void bubble_sort_preserve_indexs(double* values, int* indexes, int array_length) {
   int i, j;

   for (i = 0; i < array_length-1; i++)      
       // Last i elements are already in place   
       for (j = 0; j < array_length-i-1; j++) 
           if (values[j] > values[j+1]) {
              swap_doubles(&values[j], &values[j+1]);
              swap_ints(&indexes[j], &indexes[j+1]);
           }
}

/*Recieve desired array length n, allocate memory and return an indexes-array of same length*/
int* get_initial_indexes_array (int n) {
    int* indexes_array;
    int i;
    /*allocate memory*/
    indexes_array = (int*)malloc(n * sizeof(int));
    assert (indexes_array != NULL); //TODO: add printf off error message

    /*Insert initial indexes values to array*/
    for (i=0; i < n; i++) {
        indexes_array[i] = i;
    }

    return indexes_array;

}

/*Receives symmetric matrix and empty array, returns n*n eigenvectors matrix (including memory allocation) and fills array accordingly in-place with corresponding eigenvalues*/
Matrix* get_eigen_vectors_and_values (Matrix* originalMatrix, double* eigen_valus_array) {
    Matrix* current_A_matrix;
    Matrix* v_matrix;
    Matrix* current_P_matrix;
    double a_current_off, a_previous_off;
    int num_iterations, i, j;
    double theta,t,c,s;
    Cell* largest_non_diagonal_cell;

    largest_non_diagonal_cell = (Cell*)malloc(sizeof(Cell*));
    assert (largest_non_diagonal_cell != NULL); //TODO: add printf off error message

    /*Start with input matrix as A*/
    current_A_matrix = originalMatrix;
    /*Initialize V to be identity matrix - agnostic to matrix multiplication*/
    v_matrix = get_identity_matrix_size_n(originalMatrix -> rows);

    a_previous_off = MAXIMUM_DOUBLE;
    a_current_off = get_off_diagonal_sum_of_matrix(current_A_matrix);
    num_iterations = 0;

    while ((num_iterations < NUM_ITERATIONS) && (a_previous_off - a_current_off > EPSILON)) {
        printf("\n\nIteration number %d: \n", num_iterations+3);
        /*calculte theta,t,c,s according to formula*/
        get_cell_with_largest_value(current_A_matrix, largest_non_diagonal_cell);
        i = largest_non_diagonal_cell -> row;
        j = largest_non_diagonal_cell -> col;
        theta = (current_A_matrix -> cells[j][j] - current_A_matrix -> cells[i][i]) / (2 * current_A_matrix -> cells[i][j]);
        t = sign(theta) / (fabs(theta) + sqrt(pow(theta,2) + 1));
        c = 1 / (sqrt(pow(t,2) + 1));
        s = t * c;

        /*Get P matrix according to c & s*/
        current_P_matrix = get_rotation_matrix_for_m(current_A_matrix, largest_non_diagonal_cell, c, s);

        /*Update V according to current P*/
        v_matrix = multiply_matrices_and_free_memory(v_matrix, current_P_matrix);
        printf("Current V matrix is: \n");
        print_matrix(v_matrix);
        
        /*Preform pivot step from A to A' and update Off values*/
        current_A_matrix = preform_pivot_step_and_free_memory(current_A_matrix, i, j, c, s);
        
        a_previous_off = a_current_off;
        a_current_off = get_off_diagonal_sum_of_matrix(current_A_matrix);

        num_iterations++;
    }

    printf("\n\nFinal A matrix is:\n");
    print_matrix(current_A_matrix);

    free(largest_non_diagonal_cell); /*Largest cell data is un-needed from this point, so free the memory space allocated by get_cell_with_largest_value*/

    /*Update the eigenValues in array in-place*/
    for (i=0; i < current_A_matrix -> rows; i++) {
        eigen_valus_array[i] = current_A_matrix -> cells[i][i];
    }

    return v_matrix;
}   

/*Preform pivot from A to A' according to input A,c,s,i,j and free memory of A*/
Matrix* preform_pivot_step_and_free_memory (Matrix* A, int i, int j, double c, double s) {
    Matrix* newA;
    int r;

    newA = get_copy_of_matrix(A);

    /*iterate only on rows i,j and columns i,j */
    for (r=0; r < A -> rows; r++) {
        /*Change column i & row i*/
        newA -> cells[r][i] = (c * A -> cells[r][i]) - (s * A -> cells[r][j]);
        newA -> cells[i][r] = newA -> cells[r][i];
        /*Change column j & row j*/
        newA -> cells[r][j] = (c * A -> cells[r][j]) + (s *  A -> cells[r][i]);
        newA -> cells[j][r] = newA -> cells[r][j];
    }

    /*correct values for A'[i][i], A'[j][j] & A'[i][j]*/
    newA -> cells[i][i] = (pow(c,2) * A -> cells[i][i]) + (pow(s,2) * A -> cells[j][j]) - (2 * s * c * A -> cells[i][j]);
    newA -> cells[j][j] = (pow(s,2) * A -> cells[i][i]) + (pow(c,2) * A -> cells[j][j]) + (2 * s * c * A -> cells[i][j]);
    newA -> cells[i][j] = 0;
    newA -> cells[j][i] = 0;

    /*Free memory of original A*/
    free_matrix_memory(A);

    return newA;
}

/*Recieves symmetric matrix m and returns sum of squares of all off-diagonal elements*/
double get_off_diagonal_sum_of_matrix (Matrix* m) {
    int i,j;
    double sum = 0.0;

    /*calculte sum of squares of all elements above diagonal*/
    for (i=0; i < m -> rows; i++) {
        for (j=i+1; j < m -> rows; j++) {
            sum += pow(m -> cells[i][j], 2);
        }
    }
    /*multiply by 2 to get sum of all off-diagonal elements*/
    sum *= 2;

    return sum;
}

/*Recieves symmetric matrix m and returns rotation matrix for m*/
Matrix* get_rotation_matrix_for_m (Matrix* m, Cell* largest_non_diagonal_cell, double c, double s) {
    Matrix* rotation_matrix;
    int i,j;

    /*start with identity matrix*/
    rotation_matrix = get_identity_matrix_size_n(m -> rows);

    i = largest_non_diagonal_cell -> row;
    j = largest_non_diagonal_cell -> col;
    
    /*set correlating cells to c, s  & -s*/
    rotation_matrix -> cells[i][i] = c;
    rotation_matrix -> cells[i][j] = s;
    rotation_matrix -> cells[j][j] = c;
    rotation_matrix -> cells[j][i] = s*(-1);

    return rotation_matrix;
}