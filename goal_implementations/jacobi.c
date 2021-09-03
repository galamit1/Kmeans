//
// Created by galam on 26/07/2021.
//

#include "jacobi.h"
#include "../utils.c"
#include <assert.h>
#include <math.h>


#define EPSILON 1.0e-15

/*Recieves lnorm symmetric matrix, runs jacobi algorithm and returns n*k matrix U  containing the first k eigenvectors u1, . . . , uk of Lnorm columns. K is inferred by matrix size*/
Matrix* run_jacobi (Matrix* lnorm) {
    double* eigenValusArray;
    
}

int getEigenGapK (double* eigen_valus_array, int array_length) {
    double max_gap, current_gap;
    int max_index, current_index;

    /***TODO: solve how to preserve original values order - another copy of array? array with mapping from index to value?***/
    /*Sort the array*/
    qsort(eigen_valus_array, array_length, sizeof(double), standart_sort);
    /*Set first index to be initial max gap & index*/
    max_index = 0;
    max_gap = abs(eigen_valus_array[0] - eigen_valus_array[1]);

    /*Iterate rest of eigenValues*/
    for (current_index = 0; current_index < (array_length - 1); current_index++) {
        current_gap = abs(eigen_valus_array[current_index] - eigen_valus_array[current_index + 1]);
        /*In case of equality in the argmax of some eigengaps, use the lowest index - therefore we use > and not >=*/
        if (current_gap > max_gap) {
            max_gap = current_gap;
            max_index = current_index;
        }
    }

    /*Indexes start with 0 but we need number of eigenvalues before gap - therefore we add 1 to max_index*/
    return (max_index + 1);
}

int standart_sort (double a, double b) {
    return (a < b);
}

/*Receives symmetric matrix and empty array, returns n*n eigenvectors matrix (including memory allocation) and fills array accordingly in-place with corresponding eigenvalues*/
Matrix* getEigenVectorsAndValues (Matrix* originalMatrix, double* eigenValusArray) {
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

    a_previous_off = INFINITY*(-1);
    a_current_off = getOffDiagonalSumOfMatrix(current_A_matrix);
    num_iterations = 0;

    while ((num_iterations < 100) && (a_current_off - a_previous_off <= EPSILON)) {
        /*calculte theta,t,c,s according to formula*/
        largestNonDiagonalCell = getCellWithLargestValue(current_A_matrix);
        i = largestNonDiagonalCell -> row;
        j = largestNonDiagonalCell -> col;
        theta = (current_A_matrix -> cells[j][j] - current_A_matrix -> cells[i][i]) / (2 * current_A_matrix -> cells[i][j]);
        t = sign(theta) / (abs(theta) + sqrt(pow(theta,2) + 1));
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
        eigenValusArray[i] = current_A_matrix -> cells[i][i];
    }

    return v_matrix;
}   

/*Preform pivot from A to A' according to input A,c,s,i,j and free memory of A*/
Matrix* preformPivotStepAndFreeMemory (Matrix* A, int i, int j, double c, double s) {
    Matrix* newA;
    int r;

    newA = getZerosMatrixSizeN(A -> rows);

    /*iterate only on rows i,j and columns i,j */
    for (r=0; r < A -> rows; r++) {
        newA -> cells[r][i] = (c * A -> cells[r][i]) - (s * A -> cells[r][j]);
        newA -> cells[r][j] = (c * A -> cells[r][j]) + (s *  A -> cells[r][i]);
    }

    /*correct values for A'[i][i], A'[j][j] & A'[i][j]*/
    newA -> cells[i][i] = (pow(c,2) * A -> cells[i][i]) + (pow(s,2) * A -> cells[j][j]) - (2 * s * c * A -> cells[i][j]);
    newA -> cells[j][j] = (pow(s,2) * A -> cells[i][i]) + (pow(c,2) * A -> cells[j][j]) - (2 * s * c * A -> cells[i][j]);
    newA -> cells[i][j] = 0;

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