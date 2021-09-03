#include "spkmeans.h"


/*******************************/
/*** Main Function ***/
/*******************************/

int main(int argc, char **argv) {
    int k;
    char goal[7];
    char filePath[1000];

    FILE * file_ptr;
    errno_t err;
    int i,j;
    Matrix * points_matrix;

    int run = 0;

    if (argc != 4) {
        printf("Invalid Input!");
        exit(1);
    }
    sscanf(argv[1], "%d", &k);
    sscanf(argv[2], "%s", &goal);
    sscanf(argv[3], "%s", &filePath);

    if (k < 0) {
        printf("Invalid Input!");
        exit(1);
    }
    
    points_matrix = get_points_matrix(filePath);
//    numPoints = getNumPoints(file_ptr);
//    numCoordinates = getNumCoordinates(file_ptr); //TODO calculate it together
//    points_matrix = getMatrixFrom2DArray(points, numPoints, numCoordinates);

    if (k == 0) //TODO the huristic
    {
        k=0;
    }

    if (k >= points_matrix->rows) { //TODO > or >=?
        freeMatrixMemory(points_matrix);
        exit(1);
    }

    /***Run the function according do the goal***/
    if (strcmp(goal, "wam") == 0) {
        Matrix * wam_matrix = run_wam(points_matrix);
        printMatrix(wam_matrix);
        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(points_matrix);
        exit(0);
    }
    if (strcmp(goal, "ddg") == 0) {
        Matrix * wam_matrix = run_wam(points_matrix);
        printMatrix(wam_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        printMatrix(ddg_matrix);
        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(ddg_matrix);
        freeMatrixMemory(points_matrix);
        exit(0);
    }
    if (run == 0) {
        printf("Invalid input: goal is not in of the options");
    }

    /***clean all***/
    freeMatrixMemory(points_matrix);
}




//int main1(int argc, char **argv) {
//    /*argc = number of parameters (including program name), **argv = string array of parameters.*/
//
//    /***Variables declarations***/
//    int k;
//    char goal[7];
//    char filePath[100];
//    size_t read;
//    size_t len;
//    char * line;
//
//    /***enum Goal {SPK="spk", WAM="wam", DDG="ddg", LNORM="lnorm", JACOBI="jacobi"} goal_types;***/
//    int num_points;
//    int num_coordinates;
//    char single_line[1000]; /*Rami gave assumption each line's length is <= 1000*/
//    int line_num;
//    double **points;
//    double *point;
//    Cluster **clusters;
//    int i;
//    int j;
//    int run = 0;
//
//
//    /***Input validation***/
//    if (argc != 4) {
//        printf("Invalid input: number of parameters should be 3."); //TODO: print generic error msg
//        exit(0);
//    }
//
//    sscanf(argv[1], "%d", &k);
//    sscanf(argv[2], "%s", &goal);
//    sscanf(argv[3], "%s", &filePath);
//
//    if (k <= 0) {
//        printf("Invalid input: k should be > 0"); //TODO: print generic error msg
//        exit(0);
//    }
//
//    FILE * file = fopen(filePath, "r");
//    if (stdin == NULL) {
//        printf ("Couldn't open file"); //TODO: print generic error msg
//        exit(0);
//    }
//
//    /***Get data from file***/
//    fgets(single_line , 1000, stdin);
//    num_coordinates = get_num_coordinates(single_line);
//    num_points = 1; /*Already scanned first line*/
//    while (fgets(single_line, 1000, stdin) != NULL) {
//        num_points++;
//    }
//
//    rewind(stdin); /*Reset pointer to head of file*/
//
//    if (k > num_points) {
//        printf ("Invalid input: more clusters than points"); //TODO: print generic error msg
//        exit(0);
//    }
//
//    points = init_points(num_points, num_coordinates);
//    line_num = 0;
//    fgets(single_line, 1000, stdin); /*Rami gave assumption there's at least 1 point*/
//
//    do {
//        convert_line_to_point(points[line_num], single_line);
//        line_num++;
//    }
//    while (fgets(single_line, 1000, stdin) != NULL);
//
//    fclose(stdin);
//
//    /***Run the function according do the goal***/
//    if (strcmp(goal, "spk")) {
//        spk_run(points, k, num_points, num_coordinates);
//        run = 1;
//    }
//    if (run == 0) {
//        printf("Invalid input: goal is not in of the options");
//    }
//
//
//    /***Print results***/
//    for (i=0; i<k; i++) {
//        for (j=0; j<num_coordinates; j++) {
//            printf("%.4f", clusters[i]->curr_centroids[j]);
//            if (j != (num_coordinates-1)) {
//                printf(",");
//            }
//        }
//        printf("\n");
//    }
//
//    free_points_memory(points, num_points);
//    return 0;
//
//}

/*******************************/
/*** Function Implementation ***/
/*******************************/


int get_num_coordinates(char* sinlge_line) {
    /*
    Recieves pointer to char array containing the first line (=point),
    Returns number of coordinates in each point (we can assume all points have same amount of coordinates).
    */
    int count_coordinates = 1;
    char* pointer = sinlge_line;

    while ((pointer = strchr(pointer, ',')) != NULL) {
        pointer++;
        count_coordinates++;
    }

    return count_coordinates;
}


double** init_points (int num_points, int num_coordinates) {
    /*
   Recieves number of points and number of coordinates,
   Returns new 2D array of points with sufficient memory to store all points.
   */
    double** points;
    int i;

    /*allocate memory for 2D array of points*/
    points = (double**)malloc(sizeof(double*) * num_points);
    assert(points != NULL);
    for (i=0; i< num_points; i++) {
        points[i] = (double*)malloc(sizeof(double) * num_coordinates);
        assert(points[i] != NULL);
    }

    return points;
}


double* convert_line_to_point (double* point, char* line) {
    /*
   Recieves pointer to char array containing a line (=point) and pointer to initialized double-array (=point),
   Returns the same point after filling it with double values from the line, using string.h functio strtok.
   Reference: https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
   */
    char* tokens;
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


void free_points_memory (double** points, int num_points) {
    /*
    Recieves 2D array of pounts and number of points,
    Frees up memory used by all points.
    */
    int i;

    for (i=0; i<num_points; i++) {
        free(points[i]);
    }
}


