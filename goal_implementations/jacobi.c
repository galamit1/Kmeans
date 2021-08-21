//
// Created by galam on 26/07/2021.
//

#include "jacobi.h"
#include "../utils.c"
#include <assert.h>
#include <math.h>


#define EPSILON 0.001

/*TODO: implement iterative getEigenVectors function (with limit 0f 100 iterations), which should include getEigenValues function*/

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