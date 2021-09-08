#include "./utils.h"

#define LINE_SIZE 1001 /*Given assumption each line is no more than 1000 characters + \n character*/

/**********************/
/**TESTING ONLY - DON'T FORGET TO REMOVE**/
/**********************/

// int main () {
//     int num_points, num_coordinates;
//     FILE* file_ptr;
//     double ** points;
//     int i,j;
//     Matrix* points_matrix;

//     points = get_points_matrix("./example.txt");

//     file_ptr = fopen("./example.txt", "r");
//     assert(file_ptr != NULL);
//     num_points = get_num_points(file_ptr);
//     num_coordinates = get_num_coordinates(file_ptr);

//     points_matrix = getMatrixFrom2DArray(points, num_points, num_coordinates);

//     print_matrix(points_matrix);

//     return 0;
// }



/**********************/
/**MATRIX UTILS**/
/**********************/

/*Recieves 2D array of doubles and its dimentions, returns Matrix with same content, including memory allocation (for Matrix pointer only - not for Matrix contents)*/
Matrix* get_matrix_from_2D_array (double** points, int num_points, int num_coordinates) {
    Matrix* output;

    /*allocate memory*/
    output = (Matrix*)malloc(sizeof(Matrix));
    assert (output != NULL); //TODO: add printf off error message
    /*initialize parameters, including n*n matrix with memory*/
    output -> rows = num_points;
    output -> cols = num_coordinates;
    output -> cells = points;

    return output;
}

/*Receives Matrix and prints it - 4 decimal points*/
void print_matrix (Matrix* m) {
    int i,j;

    printf("Matrix dimentions are: %d rows, %d columns\n", m -> rows, m -> cols);
    for (i=0; i < m -> rows; i++) {
        for (j=0; j < m -> cols - 1; j++) {
             printf("%.4f,", m -> cells[i][j]);
         }
         printf("%.4f\n", m -> cells[i][j]);
     }
}


/*Receives Matrix and prints it as is*/
void print_full_matrix (Matrix* m) {
    int i,j;

    printf("Matrix dimentions are: %d rows, %d columns\n", m -> rows, m -> cols);
    for (i=0; i < m -> rows; i++) {
        for (j=0; j < m -> cols - 1; j++) {
            printf("%f,", m -> cells[i][j]);
        }
        printf("%f\n", m -> cells[i][j]);
    }
}

/*Recieves integer n, allocates memory and returns n*n zeros matrix, including memory allocation*/
Matrix* get_zeros_matrix_size_n (int n) {
    Matrix* zeros_matrix = malloc(sizeof(Matrix));
    assert (zeros_matrix != NULL); //TODO: add printf off error message
    /*initialize parameters, including n*n matrix with memory*/
    zeros_matrix -> rows = n;
    zeros_matrix -> cols = n;
    zeros_matrix -> cells = malloc(n * sizeof(double*));
    for (int i=0; i < n; i++) {
        zeros_matrix -> cells[i] = calloc(n, sizeof(double));
        assert (zeros_matrix -> cells[i] != NULL); //TODO: add printf off error message
    }
    return zeros_matrix;
}

Matrix* get_n_k_zero_matrix (int n, int k) {
    Matrix* nk_matrix;
    int i;

    nk_matrix = malloc(sizeof(Matrix)); //TODO: add printf off error message
    assert (nk_matrix != NULL);
    /*initialize parameters, including n*k matrix with memory*/
    nk_matrix -> rows = n;
    nk_matrix -> cols = k;
    nk_matrix -> cells = malloc(n * sizeof(double*));
    assert (nk_matrix -> cells != NULL); //TODO: add printf off error message
    for (i=0; i < n; i++) {
        nk_matrix -> cells[i] = calloc(n, sizeof(double));
        assert (nk_matrix -> cells[i] != NULL); //TODO: add printf off error message
    }
    return nk_matrix;
}

/*Recieves integer n, allocates memory and returns n*n identity matrix*/
Matrix* get_identity_matrix_size_n (int n) {
    int i;
    Matrix* identity;

    /*use get_zeros_matrix_size_n for memory allocation and matrix creation*/
    identity = get_zeros_matrix_size_n(n);

    /*set diagonal to 1*/
    for (i=0; i < n; i++) {
        identity -> cells[i][i] = 1;
    }

    return identity;
}

/*Recieves symmetric matrix and returns Cell with largest absolute value (above diagonal)*/
void get_cell_with_largest_value (Matrix* m, Cell* cell_pointer) {
    int i,j;
    double curr_absolute_value;
    Cell* largest_cell;
    /*Initialize cell*/
    largest_cell -> row = -1;
    largest_cell -> col = -1;
    largest_cell -> value = -MAXIMUM_DOUBLE;
    /*find largest absolute value and update largest_cell accordingly*/
    for (i=0; i < m -> rows; i++) {
        for (j=i+1; j < m -> cols; j++) {
            curr_absolute_value = fabs(m -> cells[i][j]);
            if (curr_absolute_value > largest_cell -> value) {
                largest_cell -> row = i;
                largest_cell -> col = j;
                largest_cell -> value = curr_absolute_value;
            }
        }
    }

    cell_pointer -> row = largest_cell -> row;
    cell_pointer -> col = largest_cell -> col;
    cell_pointer -> value = largest_cell -> value;
}

/*Recieves symmetric matrix and checks wether it's diagonal or not*/
int is_diagonal_matrix (Matrix* m) {
    int i,j;
    //check if there's a non-zero element above diagonal
    for (i=0; i < m -> rows; i++) {
        for (j=i+1; j < m -> cols; j++) {
            if (m -> cells[i][j] != 0)
                return 0;
        }
    }
    return 1;    
}

/*Recieves number and returns its sign*/
double sign (double num) {
    if (num >= 0)
        return 1.0;
    return -1.0;
}

/*Recieves two n*n symmetric matrices, returns their product and frees up memory used by those 2 matrices */
Matrix* multiply_matrices_and_free_memory (Matrix* m1, Matrix* m2) {
    Matrix* product;
    int i,j, k;

    product = get_zeros_matrix_size_n(m1 -> rows);

    /*TODO: can we improve this multilpication without making it too complicated?*/
    /*preform regular matrix multiplication - store result in product*/
    for (i=0; i < product -> rows; i++) {
        for (j=0; j < product -> rows; j++) {
            for (k=0; k < product -> rows; k++) {
                product -> cells[i][j] += (m1 -> cells[i][k])*(m2 -> cells[k][j]);
            }
        }
    }

    /*free memory of original matrices*/
    free_matrix_memory(m1);
    free_matrix_memory(m2);

    return product;
}

/*Recieves symmetric matrix and frees its memory*/
void free_matrix_memory (Matrix* m) {
    int n = m -> rows;
    int i;

    for (i=0; i < n; i++) {
            free(m -> cells[i]);
    }
    free(m -> cells);
    free(m);
}

/*Recieves matrix, allocates memory and returns copy of original matrix*/
Matrix* get_copy_of_matrix(Matrix* original) {
    Matrix* new_matrix;
    int i,j;

    new_matrix = get_n_k_zero_matrix(original -> rows, original -> cols);
    /*Update values do be identical*/
    for (i=0; i < original -> rows; i++) {
        for (j=0; j < original -> cols; j++) {
            new_matrix -> cells[i][j] = original -> cells[i][j];
        }
    }

    return new_matrix;
}

/*Recieves matrix, allocates memory and returns a normalized matrix*/
Matrix* normalize_matrix (Matrix* original) {
    Matrix* normalized_matrix;
    int i,j;
    double current_row_squares_sum, sqrt_current_row_sum;

    /*Start from copy of original matrix*/
    normalized_matrix = get_copy_of_matrix(original);
    /*Iterate rows*/
    for (i=0; i < original -> rows; i++) {
        /*Calculate row sum and the sqrt of it*/
        current_row_squares_sum = 0;
        for (j=0; j < original -> cols; j++) {
            current_row_squares_sum += pow(original -> cells[i][j], 2);
        }
        sqrt_current_row_sum = sqrt(current_row_squares_sum);

        /*Set normalized value for each row element accordingly*/
        for (j=0; j < original -> cols; j++) {
            normalized_matrix -> cells[i][j] = (original -> cells[i][j]) / sqrt_current_row_sum;
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

    while ((ch = fgetc(file_ptr)) != EOF) {
        if (ch == '\n') {
            count_points++;
            curr_location = ftell(file_ptr);

        }
    }
    if (curr_location == ftell(file_ptr)) { // there is an extra \n at the end of the file
        count_points --;
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
void get_points_from_file (int num_points, int num_coordinates, FILE *file_ptr, double** points) {
    int current_line;
    char single_line[LINE_SIZE];

    current_line = 0;
    /*verify file is non-empty - we have assumption that it's not*/
    assert(fgets(single_line, LINE_SIZE, file_ptr) != NULL);
    /*already read first line so use do-while*/
    do {
       single_line_to_point(points[current_line], single_line);
        current_line++;
    }  
    while (fgets(single_line, LINE_SIZE, file_ptr) != NULL); 
}

/*Recieves pointer to single point and single line (as string), inserts coordinates (as double) to point*/
void single_line_to_point (double* point, char* single_line) {
    char* token;
    char* pointer;
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
Matrix* get_points_matrix (char *path) {
    Matrix * points_matrix = malloc(sizeof(Matrix));
    FILE *file_ptr;
    double** points;
    int num_points, num_coordinates;

    if ((file_ptr = fopen(path, "r")) != NULL) {
        num_points = get_num_points(file_ptr);
        num_coordinates = get_num_coordinates(file_ptr);

        //Allocate memory for 2D-array
        points = (double**)malloc(num_points * sizeof(double*));
        if (points == NULL) {
            free(points_matrix);
            return NULL;
        }
        for (int i=0; i < num_points; i++) {
            points[i] = (double*)malloc(num_coordinates * sizeof(double));
            if (points[i] == NULL) {
                for (int j=0; j < i; j++) {
                    free(points[i]);
                }
                free(points);
                free(points_matrix);
                return NULL;
            }
        }

    /*Fill array with parsed values from files*/
    get_points_from_file(num_points , num_coordinates , file_ptr , points);
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