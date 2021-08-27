#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_SIZE 1000 /*Given assumption each line is no more than 1000 characters*/

/**********************/
/**TESTING ONLY - DON'T FORGET TO REMOVE**/
/**********************/

int main () {
    int numPoints, numCoordinates;
    FILE* fptr;

    fptr = fopen("example.txt", "r");
    numPoints = getNumPoints(fptr);
    numCoordinates = getNumCoordinates(fptr);

    printf("Number of points is: %d, Number of coordinates is: %d", numPoints, numCoordinates);

    return 0;
}


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

/*Recieves file pointer and returns number of points == number of lines in file*/
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

/*Recieves file pointer and returns number of coordinates in each line*/
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

/*Recieves number of points & coordinates along with FILE pointer and empty 2D points array, fills array with points from file*/
void getPointsFromFile (int numPoints, int numCoordinates, FILE *fptr, double** points) {
    int currentLine;
    char singleLine[LINE_SIZE];
    char* token;

    currentLine = 0;
    /*verify file is non-empty - we have assumption that it's not*/
    assert(fgets(singleLine, LINE_SIZE, fptr) != NULL);
    /*already read first line so use do-while*/
    do {
       singleLineToPoint(points[currentLine], singleLine);
        currentLine++;
    }  
    while (fgets(singleLine, LINE_SIZE, fptr) != NULL); 
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

/*Receives path to file and returns 2D matrix with all points, including memory allocation*/
double** get_points (char *path) {
    FILE *fptr;
    double** points;
    int numPoints, numCoordinates;
    int i;

    fptr = fopen(path,"r");
    if (fptr != NULL) {
        numPoints = getNumPoints(fptr);
        numCoordinates = getNumCoordinates(fptr);

        //Allocate memory for 2D-array
        points = (double**)malloc(numPoints * sizeof(double*));
        for (i=0; i < numPoints; i++) {
            points[i] = (double*)malloc(numCoordinates * sizeof(double));
            assert(points[i] != NULL); //TODO: add printf off error message
        }

    /*Fill array with parsed values from files*/
    getPointsFromFile(numPoints , numCoordinates , fptr , points);
    }

    fclose(fptr);
    
    return points;
}