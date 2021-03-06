#include "spkmeans.h"

/**********************/
/**MAIN FUNCTION**/
/**********************/

int main(int argc, char **argv) {
    /*argc = number of parameters (including program name), **argv = string array of parameters.*/
    int k;
    char goal[7];
    char filePath[1000];
    Matrix *points_matrix;

    if (argc != 4) {
        printf("Invalid Input!");
        return 1;
    }
    /* assuming the input is valid */
    sscanf(argv[1], "%d", &k);
    sscanf(argv[2], "%s", goal);
    sscanf(argv[3], "%s", filePath);

    if (k < 0) {
        printf("Invalid Input!");
        return 1;
    }

    if ((points_matrix = get_points_matrix(filePath)) == NULL) {
        printf("An Error Has Occured");
        return 1;
    }

    if (k >= points_matrix->rows) { /*so that k < number of points*/
        printf("Invalid Input!");
        free_matrix_memory(points_matrix);
        return 1;
    }
    run_functions_according_to_goal(goal, points_matrix, k);
    return 0;
}


void run_functions_according_to_goal(char *goal, Matrix *points_matrix, int k) {
    /* check what the goal is and run the right function */
    int run = 0;
    Matrix *wam_matrix;
    Matrix *ddg_matrix;
    Matrix *lnorm_matrix;
    Matrix *U_matrix;
    Matrix *T_matrix;


    if (strcmp(goal, "wam") == 0) {
        run = 1;
        wam_matrix = run_wam(points_matrix);
        print_matrix(wam_matrix);

        free_matrix_memory(wam_matrix);
    }
    if (strcmp(goal, "ddg") == 0) {
        run = 1;
        wam_matrix = run_wam(points_matrix);
        ddg_matrix = run_ddg(wam_matrix);
        print_matrix(ddg_matrix);

        free_matrix_memory(wam_matrix);
        free_matrix_memory(ddg_matrix);
    }
    if (strcmp(goal, "lnorm") == 0) {
        run = 1;
        wam_matrix = run_wam(points_matrix);
        ddg_matrix = run_ddg(wam_matrix);
        convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
        lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
        print_matrix(lnorm_matrix);

        free_matrix_memory(wam_matrix);
        free_matrix_memory(ddg_matrix);
        free_matrix_memory(lnorm_matrix);
    }
    if (strcmp(goal, "jacobi") == 0) {
        run = 1;
        U_matrix = run_jacobi(points_matrix, goal, k);

        free_matrix_memory(U_matrix);
    }
    if (strcmp(goal, "spk") == 0) {
        run = 1;
        wam_matrix = run_wam(points_matrix);
        ddg_matrix = run_ddg(wam_matrix);
        convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
        lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
        U_matrix = run_jacobi(lnorm_matrix, goal, k);
        T_matrix = normalize_matrix(U_matrix);

        if (k == 0) {
            k = T_matrix->cols;
        }
        run_spk(T_matrix, k);

        free_matrix_memory(wam_matrix);
        free_matrix_memory(ddg_matrix);
        free_matrix_memory(lnorm_matrix);
        free_matrix_memory(U_matrix);
        free_matrix_memory(T_matrix);
    }

    if (run == 0) {
        printf("Invalid input");
    }

    free_matrix_memory(points_matrix);
}

/**********************/
/**MATRIX UTILS**/
/**********************/

/*Recieves 2D array of doubles and its dimentions, returns Matrix with same content, including memory allocation (for Matrix pointer only - not for Matrix contents)*/
Matrix *get_matrix_from_2D_array(double **points, int num_points, int num_coordinates) {
    Matrix *output;

    /*allocate memory*/
    output = (Matrix *) malloc(sizeof(Matrix));
    assert (output != NULL && "An Error Has Occured");
    /*initialize parameters, including n*n matrix with memory*/
    output->rows = num_points;
    output->cols = num_coordinates;
    output->cells = points;

    return output;
}

/*Receives Matrix and prints it - 4 decimal points*/
void print_matrix(Matrix *m) {
    int i, j;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols - 1; j++) {
            if (m->cells[i][j] > -0.00005 && m->cells[i][j] < 0.0) {
                printf("%.4f,", m->cells[i][j] * (-1));
            } else {
                printf("%.4f,", m->cells[i][j]);
            }
        }
        if (m->cells[i][j] > -0.00005 && m->cells[i][j] < 0.0) {
            printf("%.4f\n", m->cells[i][j] * (-1));
        } else {
            printf("%.4f\n", m->cells[i][j]);
        }
    }
}

/*Receives Matrix and prints its transpose - 4 decimal points*/
void print_matrix_transpose(Matrix *m) {
    int i, j;
    for (j = 0; j < m->cols; j++) {
        for (i = 0; i < m->rows - 1; i++) {
            if (m->cells[i][j] > -0.00005 && m->cells[i][j] < 0.0) {
                printf("%.4f,", m->cells[i][j] * (-1));
            } else {
                printf("%.4f,", m->cells[i][j]);
            }
        }
        if (m->cells[i][j] > -0.00005 && m->cells[i][j] < 0.0) {
            printf("%.4f\n", m->cells[i][j] * (-1));
        } else {
            printf("%.4f\n", m->cells[i][j]);
        }
    }
}


/*Receives Matrix and prints it as is*/
void print_full_matrix(Matrix *m) {
    int i, j;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols - 1; j++) {
            printf("%f,", m->cells[i][j]);
        }
        printf("%f\n", m->cells[i][j]);
    }
}

/*Recieves integer n, allocates memory and returns n*n zeros matrix*/
Matrix *get_zeros_matrix_size_n(int n) {
    Matrix *zeros_matrix;
    int i;

    zeros_matrix = malloc(sizeof(Matrix));
    assert (zeros_matrix != NULL && "An Error Has Occured");
    /*initialize parameters, including n*n matrix with memory*/
    zeros_matrix->rows = n;
    zeros_matrix->cols = n;
    zeros_matrix->cells = malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) {
        zeros_matrix->cells[i] = calloc(n, sizeof(double));
        assert (zeros_matrix->cells[i] != NULL && "An Error Has Occured");
    }
    return zeros_matrix;
}

/*Recieves integers n & k, allocates memory and returns n*k zeros matrix*/
Matrix *get_n_k_zero_matrix(int n, int k) {
    Matrix *nk_matrix;
    int i;

    nk_matrix = (Matrix *) malloc(sizeof(Matrix));
    assert (nk_matrix != NULL && "An Error Has Occured");
    /*initialize parameters, including n*k matrix with memory*/
    nk_matrix->rows = n;
    nk_matrix->cols = k;
    nk_matrix->cells = malloc(n * sizeof(double *));
    assert (nk_matrix->cells != NULL && "An Error Has Occured");
    for (i = 0; i < n; i++) {
        nk_matrix->cells[i] = calloc(n, sizeof(double));
        assert (nk_matrix->cells[i] != NULL && "An Error Has Occured");
    }
    return nk_matrix;
}

/*Recieves integer n, allocates memory and returns n*n identity matrix*/
Matrix *get_identity_matrix_size_n(int n) {
    int i;
    Matrix *identity;

    /*use get_zeros_matrix_size_n for memory allocation and matrix creation*/
    identity = get_zeros_matrix_size_n(n);

    /*set diagonal to 1*/
    for (i = 0; i < n; i++) {
        identity->cells[i][i] = 1;
    }

    return identity;
}

/*Recieves symmetric matrix and returns Cell with largest absolute value (above diagonal)*/
void get_cell_with_largest_value(Matrix *m, Cell *cell_pointer) {
    int i, j;
    double curr_absolute_value;
    Cell *largest_cell;
    /*Initialize cell*/
    largest_cell = (Cell *) malloc(sizeof(Cell));
    largest_cell->row = -1;
    largest_cell->col = -1;
    largest_cell->value = -MAXIMUM_DOUBLE;
    /*find largest absolute value and update largest_cell accordingly*/
    for (i = 0; i < m->rows; i++) {
        for (j = i + 1; j < m->cols; j++) {
            curr_absolute_value = fabs(m->cells[i][j]);
            if (curr_absolute_value > largest_cell->value) {
                largest_cell->row = i;
                largest_cell->col = j;
                largest_cell->value = curr_absolute_value;
            }
        }
    }

    cell_pointer->row = largest_cell->row;
    cell_pointer->col = largest_cell->col;
    cell_pointer->value = largest_cell->value;
    free(largest_cell);
}


/*Recieves number and returns its sign*/
double sign(double num) {
    if (num >= 0)
        return 1.0;
    return -1.0;
}

/*Recieves 3 n*n symmetric matrices, multiplies first 2 and stores the result in product*/
void multiply_matrices_to_existing_pointer(Matrix *m1, Matrix *m2, Matrix *product) {
    int i, j, k;

    /*preform regular matrix multiplication - store result in product*/
    for (i = 0; i < product->rows; i++) {
        for (j = 0; j < product->rows; j++) {
            product->cells[i][j] = 0;
            for (k = 0; k < product->rows; k++) {
                product->cells[i][j] += (m1->cells[i][k]) * (m2->cells[k][j]);
            }
        }
    }
}

/*Recieves symmetric matrix and frees its memory, the cells are double**   */
void free_matrix_memory(Matrix *m) {
    int i;

    if (m != NULL) {
        if (m->cells != NULL) {
            for (i = 0; i < m->rows; i++) {
                if (m->cells[i] != NULL) {
                    free(m->cells[i]);
                }
            }
            free(m->cells);
        }
        free(m);
    }
}

/*Recieves matrix, allocates memory and returns copy of original matrix*/
Matrix *get_copy_of_matrix(Matrix *original) {
    Matrix *new_matrix;
    int i, j;

    new_matrix = get_n_k_zero_matrix(original->rows, original->cols);
    /*Update values do be identical*/
    for (i = 0; i < original->rows; i++) {
        for (j = 0; j < original->cols; j++) {
            new_matrix->cells[i][j] = original->cells[i][j];
        }
    }

    return new_matrix;
}

/*Recieves matrix, allocates memory and returns a normalized matrix*/
Matrix *normalize_matrix(Matrix *original) {
    Matrix *normalized_matrix;
    int i, j;
    double current_row_squares_sum, sqrt_current_row_sum;

    /*Start from copy of original matrix*/
    normalized_matrix = get_copy_of_matrix(original);
    /*Iterate rows*/
    for (i = 0; i < original->rows; i++) {
        /*Calculate row sum and the sqrt of it*/
        current_row_squares_sum = 0;
        for (j = 0; j < original->cols; j++) {
            current_row_squares_sum += pow(original->cells[i][j], 2);
        }
        sqrt_current_row_sum = sqrt(current_row_squares_sum);

        /*Set normalized value for each row element accordingly*/
        for (j = 0; j < original->cols; j++) {
            normalized_matrix->cells[i][j] = (original->cells[i][j]) / sqrt_current_row_sum;
        }
    }

    return normalized_matrix;
}

/**********************/
/**FILE HANDLING FUNCTIONS**/
/**********************/

/*Recieves file pointer and returns number of points == number of lines in file*/
int get_num_points(FILE *file_ptr) {
    int count_points = 1;
    char ch;
    int curr_location;

    curr_location = 0;
    while ((ch = fgetc(file_ptr)) != EOF) {
        if (ch == '\n') {
            count_points++;
            curr_location = ftell(file_ptr);
        }
    }
    if (curr_location == ftell(file_ptr)) { /* there is an extra \n at the end of the file*/
        count_points--;
    }
    rewind(file_ptr);
    return count_points;
}

/*Recieves file pointer and returns number of coordinates in each line*/
int get_num_coordinates(FILE *file_ptr) {
    int count_coordinates = 1;
    char ch;
    while ((ch = fgetc(file_ptr)) != '\n') {
        if (ch == ',') {
            count_coordinates++;
        }
    }
    rewind(file_ptr);
    return count_coordinates;
}

/*Recieves number of points & coordinates along with FILE pointer and empty 2D points array, fills array with points from file*/
void get_points_from_file(FILE *file_ptr, double **points) {
    int current_line;
    char single_line[LINE_SIZE];

    current_line = 0;
    /*verify file is non-empty - we have assumption that it's not*/
    assert(fgets(single_line, LINE_SIZE, file_ptr) != NULL && "An Error Has Occured");
    /*already read first line so use do-while*/
    do {
        single_line_to_point(points[current_line], single_line);
        current_line++;
    } while (fgets(single_line, LINE_SIZE, file_ptr) != NULL);
}

/*Recieves pointer to single point and single line (as string), inserts coordinates (as double) to point*/
void single_line_to_point(double *point, char *single_line) {
    char *token;
    char *pointer;
    int index;

    index = 0;
    token = strtok(single_line, ",");
    while (token != NULL) {
        point[index] = strtod(token, &pointer);
        token = strtok(NULL, ","); /*go to next token*/
        index++;
    }
}

/*Receives path to file and returns 2D matrix with all points, including memory allocation and closing the pointer*/
Matrix *get_points_matrix(char *path) {
    Matrix *points_matrix;
    FILE *file_ptr;
    double **points;
    int num_points, num_coordinates;
    int i, j;

    points_matrix = malloc(sizeof(Matrix));

    if ((file_ptr = fopen(path, "r")) != NULL) {
        num_points = get_num_points(file_ptr);
        num_coordinates = get_num_coordinates(file_ptr);

        /*Allocate memory for 2D-array*/
        points = (double **) malloc(num_points * sizeof(double *));
        if (points == NULL) {
            free(points_matrix);
            return NULL;
        }
        for (i = 0; i < num_points; i++) {
            points[i] = (double *) malloc(num_coordinates * sizeof(double));
            if (points[i] == NULL) {
                for (j = 0; j < i; j++) {
                    free(points[i]);
                }
                free(points);
                free(points_matrix);
                return NULL;
            }
        }

        /*Fill array with parsed values from files*/
        get_points_from_file(file_ptr, points);
    } else {
        free(points_matrix);
        return NULL;
    }

    fclose(file_ptr);
    points_matrix->cells = points;
    points_matrix->rows = num_points;
    points_matrix->cols = num_coordinates;

    return points_matrix;
}

/*******************************/
/*** WAM FUNCTIONS ***/
/*******************************/

/*Recieves matrix, allocates memory and creates its W matrix*/
Matrix *run_wam(Matrix *points) {
    Matrix *wam_matrix;
    double weight;
    int i, j;

    wam_matrix = get_zeros_matrix_size_n(points->rows);
    for (i = 0; i < points->rows; i++) {
        for (j = i + 1; j < points->rows; j++) {
            weight = calculate_weight(points->cells[i], points->cells[j], points->cols);
            wam_matrix->cells[i][j] = weight;
            wam_matrix->cells[j][i] = weight;
        }
    }
    return wam_matrix;
}

/*Recieves 2 datapoints and num_coordinates, returns the weight between both points*/
double calculate_weight(double *point1, double *point2, int num_coordinates) {
    double weight;
    int i;

    weight = 0;
    for (i = 0; i < num_coordinates; i++) {
        weight += pow(point1[i] - point2[i], 2);
    }
    return exp(-1 * (sqrt(weight) / 2));
}

/*******************************/
/*** DDG FUNCTIONS ***/
/*******************************/
/*Receives W matrix, allocates memory and creates its D matrix*/
Matrix *run_ddg(Matrix *wam) {
    Matrix *ddg_matrix;
    double row_sum;
    int i, z;

    ddg_matrix = get_zeros_matrix_size_n(wam->rows);
    for (i = 0; i < wam->rows; i++) {
        row_sum = 0;
        for (z = 0; z < wam->cols; z++) {
            row_sum += wam->cells[i][z];
        }
        ddg_matrix->cells[i][i] = row_sum;
    }
    return ddg_matrix;
}

/*Receives D matrix and converts it in-place to D^-1/2*/
void convert_ddg_with_the_pow_of_minus_half(Matrix *ddg_matrix) {
    int i;
    for (i = 0; i < ddg_matrix->rows; i++) {
        ddg_matrix->cells[i][i] = pow(ddg_matrix->cells[i][i], -0.5);
    }
}

/*******************************/
/*** LNORM FUNCTIONS ***/
/*******************************/
/*Receives W & D matrixes, computes and returns their lnorm*/
Matrix *run_lnorm(Matrix *wam, Matrix *ddg) {
    Matrix *lnorm;
    double value;
    int i, j;

    lnorm = get_identity_matrix_size_n(wam->rows);
    for (i = 0; i < wam->rows; ++i) {
        for (j = i + 1; j < wam->rows; ++j) {
            value = -1 * wam->cells[i][j] * ddg->cells[i][i] * ddg->cells[j][j];
            lnorm->cells[i][j] = value;
            lnorm->cells[j][i] = value;
        }
    }
    return lnorm;
}

/*******************************/
/*** JACOBI FUNCTIONS ***/
/*******************************/
/*Recieves lnorm symmetric matrix, runs jacobi algorithm and returns n*k matrix U  containing the first k eigenvectors u1, . . . , uk of Lnorm columns. K is inferred by matrix size*/
Matrix *run_jacobi(Matrix *lnorm, char *goal, int k) {
    Matrix *eigen_vectors_matrix;
    Matrix *final_U_matrix;
    Matrix *lnorm_copy;
    double *eigen_values_array;
    int *indexes_array;
    int n;
    int i;
    int new_k;

    n = lnorm->rows;
    /*allocate memory for eigenvalues array*/
    eigen_values_array = (double *) malloc(n * sizeof(double));
    assert (eigen_values_array != NULL && "An Error Has Occured");
    /*get eigenvectors and update eigenvalues array accordingly*/
    lnorm_copy = get_copy_of_matrix(lnorm); /*work with deep copy of lnorm to avoid freeing issues*/
    eigen_vectors_matrix = get_eigen_vectors_and_values(lnorm_copy, eigen_values_array);
    if (strcmp(goal, "jacobi") == 0) {
        /*print eigenvalues*/
        for (i = 0; i < n - 1; i++) {
            printf("%.4f,", eigen_values_array[i]);
        }
        printf("%.4f\n", eigen_values_array[n - 1]);
        /*If goal=jacobi: print eigenvectors*/
        print_matrix_transpose(eigen_vectors_matrix);
    }
    /*get initial indexes array i.e. [0,1,...,n-1] */
    indexes_array = get_initial_indexes_array(n);
    /*get eigengap k and preserve original indexes during sort*/
    new_k = get_eigen_gap_k(eigen_values_array, indexes_array, n);
    if (k==0) {
        k = new_k;
    }
    /*allocate memory for U - a n*k Matrix*/
    final_U_matrix = get_n_k_zero_matrix(n, k);

    /*point to k eigenvectors that correspond to k smallest eigenvalues*/
    for (i = 0; i < k; i++) {
        get_column_by_index(final_U_matrix, eigen_vectors_matrix, i, indexes_array[i], n);
    }

    free_matrix_memory(eigen_vectors_matrix);
    free(eigen_values_array);
    free(indexes_array);

    return final_U_matrix;
}

/*Receive 2 Matrixes, column size and column indexes for both matrixes, copy column from relevant index in m2 to relevant index in m1*/
void get_column_by_index(Matrix *m1, Matrix *m2, int m1_index, int m2_index, int column_size) {
    int i;

    for (i = 0; i < column_size; i++) {
        m1->cells[i][m1_index] = m2->cells[i][m2_index];
    }
}

/*Recieve eigenvalues array & initial indexes array, calculate k by eigengap heuristic and update indexes array along with sorting values array*/
int get_eigen_gap_k(double *eigen_valus_array, int *indexes_array, int array_length) {
    double max_gap, current_gap;
    int max_index, current_index;

    /*Sort the array and preserve indexes*/
    bubble_sort_preserve_indexs(eigen_valus_array, indexes_array, array_length);

    /*Set first index to be initial max gap & index*/
    max_index = 0;
    max_gap = fabs(eigen_valus_array[0] - eigen_valus_array[1]);

    /*Iterate rest of eigenValues*/
    for (current_index = 0; current_index < array_length / 2; current_index++) {
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

/*Recieve 2 int pointers and swap them*/
void swap_ints(int *xp, int *yp) {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

/*Recieves an array of doubles & array of indexes, bubble-sorts both arrays accordingly in-place*/
void bubble_sort_preserve_indexs(double *values, int *indexes, int array_length) {
    int i, j;

    for (i = 0; i < array_length - 1; i++)
        /* Last i elements are already in place */
        for (j = 0; j < array_length - i - 1; j++)
            if (values[j] > values[j + 1]) {
                swap_doubles(&values[j], &values[j + 1]);
                swap_ints(&indexes[j], &indexes[j + 1]);
            }
}

/*Recieve desired array length n, allocate memory and return an indexes-array of same length*/
int *get_initial_indexes_array(int n) {
    int *indexes_array;
    int i;
    /*allocate memory*/
    indexes_array = (int *) malloc(n * sizeof(int));
    assert (indexes_array != NULL && "An Error Has Occured");

    /*Insert initial indexes values to array*/
    for (i = 0; i < n; i++) {
        indexes_array[i] = i;
    }

    return indexes_array;

}

/*Receives symmetric matrix and empty array, returns n*n eigenvectors matrix (including memory allocation) and fills array accordingly in-place with corresponding eigenvalues*/
Matrix *get_eigen_vectors_and_values(Matrix *originalMatrix, double *eigen_valus_array) {
    Matrix *current_A_matrix;
    Matrix *v_matrix;
    Matrix *temp_v_matrix;
    Matrix *current_P_matrix;
    double a_current_off, a_previous_off;
    int num_iterations, i, j;
    double theta, t, c, s;
    Cell *largest_non_diagonal_cell;

    largest_non_diagonal_cell = (Cell *) malloc(sizeof(Cell));
    assert (largest_non_diagonal_cell != NULL && "An Error Has Occured");
    /*Start with input matrix as A*/
    current_A_matrix = originalMatrix;
    /*Initialize V to be identity matrix - agnostic to matrix multiplication*/
    v_matrix = get_identity_matrix_size_n(originalMatrix->rows);

    a_previous_off = MAXIMUM_DOUBLE;
    a_current_off = get_off_diagonal_sum_of_matrix(current_A_matrix);
    num_iterations = 0;

    while ((num_iterations < NUM_ITERATIONS) && (a_previous_off - a_current_off > EPSILON)) {
        /*calculte theta,t,c,s according to formula*/
        get_cell_with_largest_value(current_A_matrix, largest_non_diagonal_cell);
        i = largest_non_diagonal_cell->row;
        j = largest_non_diagonal_cell->col;
        theta = (current_A_matrix->cells[j][j] - current_A_matrix->cells[i][i]) / (2 * current_A_matrix->cells[i][j]);
        t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;

        /*Get P matrix according to c & s*/
        current_P_matrix = get_rotation_matrix_for_m(current_A_matrix, largest_non_diagonal_cell, c, s);
        temp_v_matrix = get_copy_of_matrix(v_matrix);

        /*Update V according to current P*/
        multiply_matrices_to_existing_pointer(temp_v_matrix, current_P_matrix, v_matrix);

        free_matrix_memory(temp_v_matrix);
        free_matrix_memory(current_P_matrix);

        /*Preform pivot step from A to A' and update Off values*/
        current_A_matrix = preform_pivot_step_and_free_memory(current_A_matrix, i, j, c, s);

        a_previous_off = a_current_off;
        a_current_off = get_off_diagonal_sum_of_matrix(current_A_matrix);

        num_iterations++;
    }

    free(largest_non_diagonal_cell); /*Largest cell data is un-needed from this point, so free the memory space allocated by get_cell_with_largest_value*/

    /*Update the eigenValues in array in-place*/
    for (i = 0; i < current_A_matrix->rows; i++) {
        eigen_valus_array[i] = current_A_matrix->cells[i][i];
    }

    free_matrix_memory(current_A_matrix);

    return v_matrix;
}

/*Preform pivot from A to A' according to input A,c,s,i,j and free memory of A*/
Matrix *preform_pivot_step_and_free_memory(Matrix *A, int i, int j, double c, double s) {
    Matrix *newA;
    int r;

    newA = get_copy_of_matrix(A);

    /*iterate only on rows i,j and columns i,j */
    for (r = 0; r < A->rows; r++) {
        /*Change column i & row i*/
        newA->cells[r][i] = (c * A->cells[r][i]) - (s * A->cells[r][j]);
        newA->cells[i][r] = newA->cells[r][i];
        /*Change column j & row j*/
        newA->cells[r][j] = (c * A->cells[r][j]) + (s * A->cells[r][i]);
        newA->cells[j][r] = newA->cells[r][j];
    }

    /*correct values for A'[i][i], A'[j][j] & A'[i][j]*/
    newA->cells[i][i] = (pow(c, 2) * A->cells[i][i]) + (pow(s, 2) * A->cells[j][j]) - (2 * s * c * A->cells[i][j]);
    newA->cells[j][j] = (pow(s, 2) * A->cells[i][i]) + (pow(c, 2) * A->cells[j][j]) + (2 * s * c * A->cells[i][j]);
    newA->cells[i][j] = 0;
    newA->cells[j][i] = 0;

    /*Free memory of original A*/
    free_matrix_memory(A);

    return newA;
}

/*Recieves symmetric matrix m and returns sum of squares of all off-diagonal elements*/
double get_off_diagonal_sum_of_matrix(Matrix *m) {
    int i, j;
    double sum = 0.0;

    /*calculte sum of squares of all elements above diagonal*/
    for (i = 0; i < m->rows; i++) {
        for (j = i + 1; j < m->rows; j++) {
            sum += pow(m->cells[i][j], 2);
        }
    }
    /*multiply by 2 to get sum of all off-diagonal elements*/
    sum *= 2;

    return sum;
}

/*Recieves symmetric matrix m and returns rotation matrix for m*/
Matrix *get_rotation_matrix_for_m(Matrix *m, Cell *largest_non_diagonal_cell, double c, double s) {
    Matrix *rotation_matrix;
    int i, j;

    /*start with identity matrix*/
    rotation_matrix = get_identity_matrix_size_n(m->rows);

    i = largest_non_diagonal_cell->row;
    j = largest_non_diagonal_cell->col;

    /*set correlating cells to c, s  & -s*/
    rotation_matrix->cells[i][i] = c;
    rotation_matrix->cells[i][j] = s;
    rotation_matrix->cells[j][j] = c;
    rotation_matrix->cells[j][i] = s * (-1);

    return rotation_matrix;
}

/*******************************/
/*** SPK FUNCTIONS ***/
/*******************************/
/*Receives matrix and k value, preforms spectral clustering accordingly and prints the result*/
void run_spk(Matrix *points, int k) {
    /***Variables declarations***/
    double *point;
    Cluster **clusters;
    int did_update;
    int i, j;

    /***Execute k-means algorithm***/
    clusters = init_k_clusters(points, k);
    did_update = 0;

    for (i = 0; i < MAX_ITER; i++) {
        for (j = 0; j < points->rows; j++) {
            point = points->cells[j];
            find_minimum_distance(clusters, point, k, points->cols);
        }
        did_update = update_all_clusters(clusters, k, points->cols);
        if (did_update == 1) { /*No update occured - we can break*/
            break;
        }
    }
    /***Print results***/
    for (i = 0; i < k; i++) {
        for (j = 0; j < points->cols - 1; j++) {
            if (clusters[i]->curr_centroids[j] > -0.00005 && clusters[i]->curr_centroids[j] < 0.0) {
                printf("%.4f,", clusters[i]->curr_centroids[j] * (-1));
            } else {
                printf("%.4f,", clusters[i]->curr_centroids[j]);
            }
        }
        if (clusters[i]->curr_centroids[j] > -0.00005 && clusters[i]->curr_centroids[j] < 0.0) {
            printf("%.4f\n", clusters[i]->curr_centroids[j] * (-1));
        } else {
            printf("%.4f\n", clusters[i]->curr_centroids[j]);
        }
    }

    free_clusters_memory(clusters, k);
}

/*******************************/
/*** SPK & kmeans++ FUNCTIONS ***/
/*******************************/

double **init_points(int num_points, int num_coordinates) {
    /*
   Recieves number of points and number of coordinates,
   Returns new 2D array of points with sufficient memory to store all points.
   */
    double **points;
    int i;

    /*allocate memory for 2D array of points*/
    points = (double **) malloc(sizeof(double *) * num_points);
    assert(points != NULL && "An Error Has Occured");
    for (i = 0; i < num_points; i++) {
        points[i] = (double *) malloc(sizeof(double) * num_coordinates);
        assert(points[i] != NULL && "An Error Has Occured");
    }

    return points;
}


double *convert_line_to_point(double *point, char *line) {
    /*
   Recieves pointer to char array containing a line (=point) and pointer to initialized double-array (=point),
   Returns the same point after filling it with double values from the line, using string.h functio strtok.
   Reference: https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
   */
    char *tokens;
    int index;

    index = 0;

    tokens = strtok(line, ",");

    while (tokens != NULL) {
        point[index] = atof(tokens);
        tokens = strtok(NULL, ",");
        index++;
    }

    return point;
}


Cluster *make_cluster(const double *point, const int num_coordinates, const int index) {
    /*
    Recieves a point, number of coordinates and a specific cluster index,
    Returns new Cluster with sufficient memory space to store both the current centroid and the previous centroid.
    */
    int i;
    Cluster *cluster;

    /*allocate memory*/
    cluster = (Cluster *) malloc(sizeof(Cluster));
    assert (cluster != NULL && "An Error Has Occured");
    /*initialize parameters*/
    cluster->cluster_size = 0;
    cluster->cluster_index = index;
    /*initialize centroids*/
    cluster->curr_centroids = (double *) malloc(sizeof(double) * num_coordinates);
    assert (cluster->curr_centroids != NULL && "An Error Has Occured");
    cluster->prev_centroids = (double *) malloc(sizeof(double) * num_coordinates);
    assert (cluster->prev_centroids != NULL && "An Error Has Occured");
    cluster->sum = (double *) malloc(sizeof(double) * num_coordinates);
    assert (cluster->sum != NULL && "An Error Has Occured");

    for (i = 0; i < num_coordinates; i++) {
        cluster->curr_centroids[i] = point[i];
        cluster->sum[i] = 0;
    }

    return cluster;
}

void add_point_to_cluster(Cluster *cluster, const double *point, int num_coordinates) {
    /*
    Recieves pointer to a Cluster, a point and number of coordinates,
    Adds the given point to the given cluster.
    */
    int i;

    cluster->cluster_size += 1;
    for (i = 0; i < num_coordinates; i++) {
        cluster->sum[i] += point[i];
    }
}

void find_minimum_distance(Cluster **clusters, double *point, int k, int num_coordinates) {
    /*
    Recieves 2D Clsuter array, pointer to a point, k and number of coordinates,
    Checks which centroid is closest to the given point,
    Adds point to closest cluster.
    */
    int i;
    int index_min;
    double distance_min;
    double current_distance;

    /*initialize distance_min to largest double value ~ INF*/
    distance_min = 1.7976931348623155e+308;
    index_min = -1;

    for (i = 0; i < k; i++) {
        current_distance = get_distance(clusters[i], point, num_coordinates);
        if (current_distance < distance_min) {
            distance_min = current_distance;
            index_min = clusters[i]->cluster_index;
        }
    }

    add_point_to_cluster(clusters[index_min], point, num_coordinates);
}

double update_cluster_centroids_and_sum(Cluster *cluster, int num_coordinates) {
    /*
    Recieves pointer to a changed cluster and number of coordinates,
    Updates cluster accordingly,
    Returns the distance between previous centroid and current centroid.
    */
    int i;
    double distance;

    for (i = 0; i < num_coordinates; i++) {
        cluster->prev_centroids[i] = cluster->curr_centroids[i];
        cluster->curr_centroids[i] = cluster->sum[i] / cluster->cluster_size;
        cluster->sum[i] = 0;
    }

    cluster->cluster_size = 0;
    distance = get_distance(cluster, cluster->prev_centroids, num_coordinates);
    return distance;
}

int update_all_clusters(Cluster **clusters, int k, int num_coordinates) {
    /*
    Recieves 2D array of clusters, k and number of coordinates,
    Upadates all clusters using update_cluster_centroids_and_sum while counting total change in distance,
    Returns 1 if no update was made, 0 otherwise.
    */
    int i;
    double total_update;

    total_update = 0.0;

    for (i = 0; i < k; i++) {
        total_update += update_cluster_centroids_and_sum(clusters[i], num_coordinates);
    }

    if (total_update == 0.0) {
        /*no update was made*/
        return 1;
    }
    return 0;
}

void free_clusters_memory(Cluster **clusters, int k) {
    /*
    Recieves 2D array of clusters and k,
    Frees up memory used by all clusters.
    */
    int i;
    Cluster *curr_cluster;

    for (i = 0; i < k; i++) {
        curr_cluster = clusters[i];
        if (curr_cluster != NULL) {
            free(curr_cluster->curr_centroids);
            free(curr_cluster->prev_centroids);
            free(curr_cluster->sum);
            free(curr_cluster);
        }
    }
    if (clusters != NULL) {
        free(clusters);
    }
}

void free_points_memory(double **points, int num_points) {
    /*
    Recieves 2D array of pounts and number of points,
    Frees up memory used by all points.
    */
    int i;

    for (i = 0; i < num_points; i++) {
        free(points[i]);
    }
}


double get_distance(Cluster *cluster, const double *point, int num_coordinates) {
    /*
    Recieves pointer to cluster, pointer to point and number of coordinates,
    Returns euclidean distance of point from cluster.
    */
    double distance;
    double to_add;
    int i;

    distance = 0;
    for (i = 0; i < num_coordinates; i++) {
        to_add = (cluster->curr_centroids[i] - point[i]);
        distance += (to_add * to_add);
    }
    return distance;
}

/*Recieves pointer to 2D array of points, k and number of coordinates,Returns new 2D array of Clusters with sufficient memory, initialized with the first k points.*/
Cluster **init_k_clusters(Matrix *points, int k) {
    Cluster **clusters;
    int i;

    clusters = (Cluster **) malloc(sizeof(Cluster *) * k);
    assert (clusters != NULL && "An Error Has Occured");

    for (i = 0; i < k; i++) {
        clusters[i] = make_cluster(points->cells[i], points->cols, i);
    }

    return clusters;
}

/* usage in spkmeansmodule */
Cluster **python_init_k_clusters(int k) {
    /*
    Recieves pointer to 2D array of points, k, number of coordinates and 2D array of initial indexes,
    Returns new 2D array of Clusters with sufficient memory, initialized with the first k points.
    */
    Cluster **clusters;

    clusters = malloc(sizeof(Cluster *) * k);
    assert (clusters != NULL);

    return clusters;
}