//
// Created by galam on 26/07/2021.
//

#include "jacobi.h"
#include <assert.h>
#include <math.h>


#define EPSILON 0.001

typedef struct Matrix {
    double** cells;
    int order;
    //add rows & cols instead of order
} Matrix ;

typedef struct Cell {
    int row;
    int col;
    double value;
} Cell;

/*TODO: implement iterative getEigenVectors function (with limit 0f 100 iterations), which should include getEigenValues function*/
/*TODO: utils functions to separate file*/

/*Preform pivot from A to A' according to input A,c,s,i,j and free memory of A*/
Matrix* preformPivotStepAndFreeMemory (Matrix* A, int i, int j, double c, double s) {
    Matrix* newA;
    int r;

    newA = getZerosMatrixSizeN(A -> order);

    /*iterate only on rows i,j and columns i,j */
    for (r=0; r < A -> order; r++) {
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
    for (i=0; i < m -> order; i++) {
        for (j=i+1; j < m -> order; j++) {
            sum += pow(m -> cells[i][j], 2);
        }
    }
    /*multiply by 2 to get sum of all off-diagonal elements*/
    sum *= 2;

    return sum;

}

/*Recieves symmetric matrix m and returns rotation matrix for m*/
Matrix* getRotationMatrixForM (Matrix* m) {
    Matrix* rotationMatrix;
    Cell* largestNonDiagonalCell;
    int i,j;
    double theta,t,c,s;

    /*start with identity matrix*/
    rotationMatrix = getIdentitiyMatrixSizeN(m -> order);
    
    /*calculte theta,t,c,s according to formula*/
    largestNonDiagonalCell = getIdentitiyMatrixSizeN(m);
    i = largestNonDiagonalCell -> row;
    j = largestNonDiagonalCell -> col;
    theta = (m -> cells[j][j] - m -> cells[i][i]) / (2 * m -> cells[i][j]);
    t = sign(theta) / (abs(theta) + sqrt(pow(theta,2) + 1));
    c = 1 / (sqrt(pow(t,2) + 1));
    s = t * c;
    
    /*set correlating cells to c, s  & -s*/
    rotationMatrix -> cells[i][i] = c;
    rotationMatrix -> cells[i][j] = s;
    rotationMatrix -> cells[j][j] = c;
    rotationMatrix -> cells[j][i] = s*(-1);

    return rotationMatrix;
}

/*Allocates memory and returns n*n zeros matrix*/
Matrix* getZerosMatrixSizeN (int n) {
    int i;
    Matrix* zerosMatrix;

    /*allocate memory*/
    zerosMatrix = malloc(sizeof(Matrix));
    assert (zerosMatrix != NULL); //TODO: add printf off error message
    /*initialize parameters, including n*n matrix with memory*/
    zerosMatrix -> order = n;
    zerosMatrix -> cells = calloc(n, sizeof(double*));
    for (i=0; i < n; i++) {
        zerosMatrix -> cells[i] = calloc(n, sizeof(double));
        assert (zerosMatrix -> cells[i] != NULL); //TODO: add printf off error message
    }

    return zerosMatrix;
}

/*Allocates mempry and returns n*n identity matrix*/
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
Cell* getCellWithLargestValue (Matrix* m) {
    int i,j, currAbsoluteValue;
    Cell* largestCell;
    
    /*Initialize cell*/
    largestCell = malloc(sizeof(Cell));
    assert (largestCell != NULL); //TODO: add printf off error message
    largestCell -> row = -1;
    largestCell -> col = -1;
    largestCell -> value = -INFINITY;

    /*find largest absolute value and update largestCell accordingly*/
    for (i=0; i < m -> order; i++) {
        for (j=i+1; j < m -> order; j++) {
            currAbsoluteValue = abs(m -> cells[i][j]);
            if (currAbsoluteValue > largestCell -> value) {
                largestCell -> row = i;
                largestCell -> col = j;
                largestCell -> value = currAbsoluteValue;
            }
        }
    }

    return largestCell;
}

/*Recieves symmetric matrix and checks wether it's diagonal or not*/
int isDiagonalMatrix (Matrix* m) {
    int i,j;
    //check if there's a non-zero element above diagonal
    for (i=0; i < m -> order; i++) {
        for (j=i+1; j < m -> order; j++) {
            if (m -> cells[i][j] != 0)
                return 0;
        }
    }
    return 1;    
}

/*Recieves number and returns its sign*/
int sign (double num) {
    if (num >= 0)
        return 1;
    return -1;
}

/*Recieves two n*n symmetric matrices, returns their product and frees up memory used by those 2 matrices */
Matrix* multiplyMatricesAndFreeMemory (Matrix* m1, Matrix* m2) {
    Matrix* product;
    int i,j, k;

    product = getZerosMatrixSizeN(m1 -> order);

    /*TODO: can we improve this multilpication without making it too complicated?*/
    /*preform regular matrix multiplication - store result in product*/
    for (i=0; i < product -> order; i++) {
        for (j=0; j < product -> order; j++) {
            for (k=0; k < product -> order; k++) {
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
    int n = m -> order;
    int i;

    for (i=0; i < n; i++) {
            free(m -> cells[i]);
    }
}
