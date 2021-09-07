#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"


/*******************************/
/*** Main Function ***/
/*******************************/
int run_functions_according_to_goal(char * goal, Matrix * points_matrix) {
    if (strcmp(goal, "wam") == 0) {
        Matrix * wam_matrix = run_wam(points_matrix);
        printMatrix(wam_matrix);

        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(points_matrix);
        exit(0);
    }
    if (strcmp(goal, "ddg") == 0) {
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        printMatrix(ddg_matrix);

        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(ddg_matrix);
        freeMatrixMemory(points_matrix);
        exit(0);
    }
    if (strcmp(goal, "lnorm") == 0) {
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
        Matrix * lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
        printFullMatrix(lnorm_matrix);

        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(ddg_matrix);
        freeMatrixMemory(lnorm_matrix);
        freeMatrixMemory(points_matrix);
        exit(0);
    }
    if (strcmp(goal, "jacobi") == 0) {
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
        Matrix * lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
        Matrix * jacobi_matrix = run_jacobi(lnorm_matrix);
        printMatrix(jacobi_matrix);

        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(ddg_matrix);
        freeMatrixMemory(lnorm_matrix);
        freeMatrixMemory(jacobi_matrix);
        freeMatrixMemory(points_matrix);
        exit(0);

    }
    if (run == 0) {
        printf("Invalid input: goal is not in of the options");
    }

    /***clean all***/
    freeMatrixMemory(points_matrix);
}



int main(int argc, char **argv) {
    int k;
    char goal[7];
    char filePath[1000];

    Matrix * points_matrix;

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

    if ((points_matrix = get_points_matrix(filePath)) == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }

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
    run_functions_according_to_goal(goal, points_matrix, k);
}


void run_functions_according_to_goal(char * goal, Matrix * points_matrix, int k) {
    int run = 0;

    if (strcmp(goal, "wam") == 0) {
        run = 1;
        Matrix * wam_matrix = run_wam(points_matrix);
        printMatrix(wam_matrix);

        freeMatrixMemory(wam_matrix);
    }
    if (strcmp(goal, "ddg") == 0) {
        run = 1;
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        printMatrix(ddg_matrix);

        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(ddg_matrix);
    }
    if (strcmp(goal, "lnorm") == 0) {
        run = 1;
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
        Matrix * lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
        printFullMatrix(lnorm_matrix);

        freeMatrixMemory(wam_matrix);
        freeMatrixMemory(ddg_matrix);
        freeMatrixMemory(lnorm_matrix);
    }
    if (strcmp(goal, "spk") == 0) {
        run = 1;
        run_spk(points_matrix, k);
    }

    if (run == 0) {
        printf("Invalid input");
    }

    /***clean all***/
    freeMatrixMemory(points_matrix);
}


static PyObject* c_spkmeans(PyObject *self, PyObject *args) {
    /*Define variables to receive from user*/
    int k;
    int max_iter = 300;
    char *goal;
    int num_points;
    int num_coordinates;
    PyObject *data_points;
    double **points;
    Matrix * points_matrix;


    /* Parse arguments from Python */
    if ((!PyArg_ParseTuple(args, "Osii", &data_points, &goal, &num_points, &num_coordinates))) {
        printf("An Error Has Occured");
        return NULL; /*In the CPython API, Null is never a valid value for a PyObject* - so it signals an error*/
    }

    /*Verify that data_points & initial_indexes are python lists*/
    if (!PyList_Check(data_points)) {
        printf("An Error Has Occured");
        return NULL;
    }

    /*Load points from PyObject into C format of 2D-array points*/
    points = init_points(num_points, num_coordinates);

    if (python_list_of_lists_to_2D_array(data_points, points) != 0) {
        return NULL;
    }
    if ((points_matrix = malloc(sizeof(Matrix))) == NULL) {
        printf("An Error Has Occured");
    }

    points_matrix->rows = num_points;
    points_matrix->cols = num_coordinates;
    points_matrix->cells = points;
    run_functions_according_to_goal(goal, points_matrix, k);
}

/*******************************/
/*** Module setup ***/
/*******************************/

static PyMethodDef _methods[] = {
        {"fit",
                (PyCFunction) c_spkmeans,
                     METH_VARARGS,
                PyDoc_STR("run functions according to goal"),
        },
        {NULL, NULL, 0, NULL} /*Sentinel*/
};

static struct PyModuleDef _moduledef  = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans",
        NULL,
        -1,
        _methods
};

PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}


