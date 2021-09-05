//
// Created by galam on 26/07/2021.
//

#include "./jacobi.h"

#define EPSILON 1.0e-15
#define NUM_ITERATIONS 100

/**********************/
/**TESTING ONLY - DON'T FORGET TO REMOVE**/
/**********************/

 int main () {
     Matrix* lnorm;
     Matrix* lnorm_copy;
     Matrix* jacobi;

     lnorm = get_points_matrix("./lnorm-example.txt");
    //  printf("Original: \n");
    //  printMatrix(lnorm); 

    jacobi = run_jacobi(lnorm);
    printMatrix(jacobi);
    
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
    int i;

    n = lnorm -> rows;
    /*allocate memory for eigenvalues array*/
    eigen_valus_array = (double*)malloc(n * sizeof(double));
    assert (eigen_valus_array != NULL); //TODO: add printf off error message
    /*get eigenvectors and update eigenvalues array accordingly*/
    eigen_vectors_matrix = getEigenVectorsAndValues(lnorm, eigen_valus_array);
    /*******DELETE******/
    printf("Eigenvectors matrix is: \n");
    printMatrix(eigen_vectors_matrix);
    /*******DELETE******/
    /*get initial indexes array i.e. [0,1,...,n-1] */
    indexes_array = get_initial_indexes_array(n);
    /*get eigengap k and preserve original indexes during sort*/
    k = getEigenGapK(eigen_valus_array, indexes_array, n);

    /*allocate memory for U - a n*k Matrix*/
    final_U_matrix = malloc(sizeof(Matrix));
    assert (final_U_matrix != NULL); //TODO: add printf off error message
    /*initialize parameters, including n*k matrix with memory*/
    final_U_matrix -> rows = n;
    final_U_matrix -> cols = k;
    final_U_matrix -> cells = malloc(n * sizeof(double*));
    assert (final_U_matrix -> cells != NULL); //TODO: add printf off error message
    for (i=0; i < n; i++) {
        final_U_matrix -> cells[i] = (double*)malloc(k * sizeof(double*));
        assert (final_U_matrix -> cells[i] != NULL); //TODO: add printf off error message
    }

    /*point to k eigenvectors that correspond to k smallest eigenvalues*/
    for (i=0; i < k; i++) {
        get_column_by_index(final_U_matrix, eigen_vectors_matrix, indexes_array[i], i, n);
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
int getEigenGapK (double* eigen_valus_array, int* indexes_array, int array_length) {
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
Matrix* getEigenVectorsAndValues (Matrix* originalMatrix, double* eigen_valus_array) {
    Matrix* current_A_matrix;
    Matrix* v_matrix;
    Matrix* current_P_matrix;
    double a_current_off, a_previous_off;
    int num_iterations, i, j;
    double theta,t,c,s;
    Cell* largestNonDiagonalCell;

    /*Start with input matrix as A*/
    current_A_matrix = originalMatrix;
    /*Initialize V to be identity matrix - agnostic to matrix multiplication*/
    v_matrix = getIdentitiyMatrixSizeN(originalMatrix -> rows);

    a_previous_off = MAXIMUM_DOUBLE;
    a_current_off = getOffDiagonalSumOfMatrix(current_A_matrix);
    num_iterations = 0;

    while ((num_iterations < NUM_ITERATIONS) && (a_previous_off - a_current_off > EPSILON)) {
        /*calculte theta,t,c,s according to formula*/
        largestNonDiagonalCell = getCellWithLargestValue(current_A_matrix);
        i = largestNonDiagonalCell -> row;
        j = largestNonDiagonalCell -> col;
        theta = (current_A_matrix -> cells[j][j] - current_A_matrix -> cells[i][i]) / (2 * current_A_matrix -> cells[i][j]);
        t = sign(theta) / (fabs(theta) + sqrt(pow(theta,2) + 1));
        c = 1 / (sqrt(pow(t,2) + 1));
        s = t * c;

        /*Get P matrix according to c & s*/
        current_P_matrix = getRotationMatrixForM(current_A_matrix, largestNonDiagonalCell, c, s);

        /*Update V according to current P*/
        v_matrix = multiplyMatricesAndFreeMemory(v_matrix, current_P_matrix);
        
        /*Preform pivot step from A to A' and update Off values*/
        current_A_matrix = preformPivotStepAndFreeMemory(current_A_matrix, i, j, c, s);
        
        a_previous_off = a_current_off;
        a_current_off = getOffDiagonalSumOfMatrix(current_A_matrix);

        num_iterations++;
    }

    /*Update the eigenValues in array in-place*/
    for (i=0; i < current_A_matrix -> rows; i++) {
        eigen_valus_array[i] = current_A_matrix -> cells[i][i];
    }

    return v_matrix;
}   

/*Preform pivot from A to A' according to input A,c,s,i,j and free memory of A*/
Matrix* preformPivotStepAndFreeMemory (Matrix* A, int i, int j, double c, double s) {
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
    freeMatrixMemory(A);

    return newA;
}

/*Recieves symmetric matrix m and returns sum of squares of all off-diagonal elements*/
double getOffDiagonalSumOfMatrix (Matrix* m) {
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
Matrix* getRotationMatrixForM (Matrix* m, Cell* largestNonDiagonalCell, double c, double s) {
    Matrix* rotationMatrix;
    int i,j;

    /*start with identity matrix*/
    rotationMatrix = getIdentitiyMatrixSizeN(m -> rows);

    i = largestNonDiagonalCell -> row;
    j = largestNonDiagonalCell -> col;
    
    /*set correlating cells to c, s  & -s*/
    rotationMatrix -> cells[i][i] = c;
    rotationMatrix -> cells[i][j] = s;
    rotationMatrix -> cells[j][j] = c;
    rotationMatrix -> cells[j][i] = s*(-1);

    return rotationMatrix;
}