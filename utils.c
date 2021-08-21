#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**********************/
/**STRUCT DEFINITIONS**/
/**********************/

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

/**********************/
/**MATRIX UTILS**/
/**********************/

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

/**********************/
/**FILE HANDLING FUNCTIONS**/
/**********************/

int getNumPoints(FILE *fptr) {
    int countPoints = 1;
    char ch;
    while (ch = fgetc(fptr) != EOF) {
        if (ch == '\n') {
            countPoints++;
        }
    }
    rewind(fptr);
    return countPoints;
}

int getNumCoordinates(FILE *fptr) {
    int countCoordinates = 1;
    char ch;
    while (ch = fgetc(fptr) != '\n') {
        if (ch == ',') {
            countCoordinates++;
        }
    }
    rewind(fptr);
    return countCoordinates;
}

/*TODO*/
void getPointsFromFile (int numPoints, int numCoordinates, FILE *fptr, double** points) {
    for (int i=0; i<numPoints; i++) {
        fscanf(fptr , "%lf,%lf,%lf\n", points[i * numPoints] , points[i * numPoints + 1] , points[i * numPoints + 2]); //TODO: read points from line dynamically, Reference: https://stackoverflow.com/questions/43286609/how-to-use-fscanf-to-read-doubles-in-c 
    }
    
}

double** get_points (char *filename) {
    FILE *fptr;
    double** points;
    int numPoints, numCoordinates;
    int i;

    fptr = fopen(filename,"r");
    if (fptr != NULL) {
        numPoints = getNumPoints(fptr);
        numCoordinates = getNumCoordinates(fptr);

        //Allocate memory for 2D-array
        points = malloc(numPoints * sizeof(double*));
        for (i=0; i < numPoints; i++) {
            points[i] = malloc(numCoordinates * sizeof(double));
            assert(points[i] != NULL); //TODO: add printf off error message
        }

    
    //Fill array with parsed values from files
    getPointsFromFile(numPoints , numCoordinates , fptr , points);

    }

    fclose(fptr);
    
    return points;
}

