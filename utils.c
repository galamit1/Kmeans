#include "./utils.h"

#define LINE_SIZE 1001 /*Given assumption each line is no more than 1000 characters + \n character*/

/**********************/
/**TESTING ONLY - DON'T FORGET TO REMOVE**/
/**********************/

// int main () {
//     int numPoints, numCoordinates;
//     FILE* file_ptr;
//     double ** points;
//     int i,j;
//     Matrix* points_matrix;

//     points = get_points_matrix("./example.txt");

//     file_ptr = fopen("./example.txt", "r");
//     assert(file_ptr != NULL);
//     numPoints = getNumPoints(file_ptr);
//     numCoordinates = getNumCoordinates(file_ptr);

//     points_matrix = getMatrixFrom2DArray(points, numPoints, numCoordinates);

//     printMatrix(points_matrix);

//     return 0;
// }



/**********************/
/**MATRIX UTILS**/
/**********************/

/*Recieves 2D array of doubles and its dimentions, returns Matrix with same content, including memory allocation (for Matrix pointer only - not for Matrix contents)*/
Matrix* getMatrixFrom2DArray (double** points, int numPoints, int numCoordinates) {
    Matrix* output;
    int i;

    /*allocate memory*/
    output = (Matrix*)malloc(sizeof(Matrix));
    assert (output != NULL); //TODO: add printf off error message
    /*initialize parameters, including n*n matrix with memory*/
    output -> rows = numPoints;
    output -> cols = numCoordinates;
    output -> cells = points;

    return output;
}

/*Receives Matrix and prints it - 4 decimal points*/
void printMatrix (Matrix* m) {
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
void printFullMatrix (Matrix* m) {
    double** cells;
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
Matrix* getZerosMatrixSizeN (int n) {
    Matrix* zerosMatrix = malloc(sizeof(Matrix));
    assert (zerosMatrix != NULL); //TODO: add printf off error message
    /*initialize parameters, including n*n matrix with memory*/
    zerosMatrix -> rows = n;
    zerosMatrix -> cols = n;
    zerosMatrix -> cells = malloc(n * sizeof(double*));
    for (int i=0; i < n; i++) {
        zerosMatrix -> cells[i] = calloc(n, sizeof(double));
        assert (zerosMatrix -> cells[i] != NULL); //TODO: add printf off error message
    }
    return zerosMatrix;
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
Matrix* getIdentitiyMatrixSizeN (int n) {
    int i;
    Matrix* identity;

    /*use getZerosMatrixSizeN for memory allocation and matrix creation*/
    identity = getZerosMatrixSizeN(n);

    /*set diagonal to 1*/
    for (i=0; i < n; i++) {
        identity -> cells[i][i] = 1;
    }

    return identity;
}

/*Recieves symmetric matrix and returns Cell with largest absolute value (above diagonal)*/
void getCellWithLargestValue (Matrix* m, Cell* cell_pointer) {
    int i,j;
    double currAbsoluteValue;
    Cell* largestCell;
    /*Initialize cell*/
    largestCell -> row = -1;
    largestCell -> col = -1;
    largestCell -> value = -MAXIMUM_DOUBLE;
    /*find largest absolute value and update largestCell accordingly*/
    for (i=0; i < m -> rows; i++) {
        for (j=i+1; j < m -> cols; j++) {
            currAbsoluteValue = fabs(m -> cells[i][j]);
            if (currAbsoluteValue > largestCell -> value) {
                largestCell -> row = i;
                largestCell -> col = j;
                largestCell -> value = currAbsoluteValue;
            }
        }
    }

    cell_pointer -> row = largestCell -> row;
    cell_pointer -> col = largestCell -> col;
    cell_pointer -> value = largestCell -> value;
}

/*Recieves symmetric matrix and checks wether it's diagonal or not*/
int isDiagonalMatrix (Matrix* m) {
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
Matrix* multiplyMatricesAndFreeMemory (Matrix* m1, Matrix* m2) {
    Matrix* product;
    int i,j, k;

    product = getZerosMatrixSizeN(m1 -> rows);

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
    freeMatrixMemory(m1);
    freeMatrixMemory(m2);

    return product;
}

/*Recieves symmetric matrix and frees its memory*/
void freeMatrixMemory (Matrix* m) {
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
    int n;
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
int getNumPoints(FILE *file_ptr) {
    int countPoints = 1;
    char ch;
    int curr_location;

    while ((ch = fgetc(file_ptr)) != EOF) {
        if (ch == '\n') {
            countPoints++;
            curr_location = ftell(file_ptr);

        }
    }
    if (curr_location == ftell(file_ptr)) { // there is an extra \n at the end of the file
        countPoints --;
    }
    rewind(file_ptr);
    return countPoints;
}

/*Recieves file pointer and returns number of coordinates in each line*/
int getNumCoordinates(FILE *file_ptr) {
    int countCoordinates = 1;
    char ch;
    while ((ch = fgetc(file_ptr)) != '\n') {
        if (ch == ',') {
            countCoordinates++;
        }
    }
    rewind(file_ptr);
    return countCoordinates;
}

/*Recieves number of points & coordinates along with FILE pointer and empty 2D points array, fills array with points from file*/
void getPointsFromFile (int numPoints, int numCoordinates, FILE *file_ptr, double** points) {
    int currentLine;
    char singleLine[LINE_SIZE];
    char* token;

    currentLine = 0;
    /*verify file is non-empty - we have assumption that it's not*/
    assert(fgets(singleLine, LINE_SIZE, file_ptr) != NULL);
    /*already read first line so use do-while*/
    do {
       singleLineToPoint(points[currentLine], singleLine);
        currentLine++;
    }  
    while (fgets(singleLine, LINE_SIZE, file_ptr) != NULL); 
}

/*Recieves pointer to single point and single line (as string), inserts coordinates (as double) to point*/
void singleLineToPoint (double* point, char* singleLine) {
    char* token;
    char* pointer;
    int index;

    index = 0;
    token = strtok(singleLine, ",");
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
    int numPoints, numCoordinates;

    if ((file_ptr = fopen(path, "r")) != NULL) {
        numPoints = getNumPoints(file_ptr);
        numCoordinates = getNumCoordinates(file_ptr);

        //Allocate memory for 2D-array
        points = (double**)malloc(numPoints * sizeof(double*));
        if (points == NULL) {
            free(points_matrix);
            return NULL;
        }
        for (int i=0; i < numPoints; i++) {
            points[i] = (double*)malloc(numCoordinates * sizeof(double));
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
    getPointsFromFile(numPoints , numCoordinates , file_ptr , points);
    } else {
        free(points_matrix);
        return NULL;
    }

    fclose(file_ptr);
    points_matrix->cells = points;
    points_matrix->rows = numPoints;
    points_matrix->cols = numCoordinates;
    
    return points_matrix;
}