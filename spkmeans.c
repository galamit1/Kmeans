#define PY_SSIZE_T_CLEAN
#include <python.h>
#include "spkmeans.h"


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

    if ((k == 0) && ((strcmp(goal, "spk")) || (strcmp(goal, "jacobi"))))
    {
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
        Matrix * lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
        Matrix* U_matrix = run_jacobi(lnorm_matrix);
        Matrix* T_matrix = normalize_matrix(U_matrix);

        if (strcmp(goal, "spk") == 0) {
            k = T_matrix->cols;
            run_spk(T_matrix, k);
        }
        if (strcmp(goal, "jacobi") == 0) {
            print_matrix(T_matrix);
        }

        free_matrix_memory(wam_matrix);
        free_matrix_memory(ddg_matrix);
        free_matrix_memory(lnorm_matrix);
        free_matrix_memory(U_matrix);
        free_matrix_memory(T_matrix);
        free_matrix_memory(points_matrix);
        exit(0);
    }

    if (k >= points_matrix->rows) { // so that k < number of points
        printf("Invalid Input!");
        free_matrix_memory(points_matrix);
        exit(1);
    }
    run_functions_according_to_goal(goal, points_matrix, k);
}


void run_functions_according_to_goal(char * goal, Matrix * points_matrix, int k) {
    int run = 0;

    if (strcmp(goal, "wam") == 0) {
        run = 1;
        Matrix * wam_matrix = run_wam(points_matrix);
        print_matrix(wam_matrix);

        free_matrix_memory(wam_matrix);
    }
    if (strcmp(goal, "ddg") == 0) {
        run = 1;
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        print_matrix(ddg_matrix);

        free_matrix_memory(wam_matrix);
        free_matrix_memory(ddg_matrix);
    }
    if (strcmp(goal, "lnorm") == 0) {
        run = 1;
        Matrix * wam_matrix = run_wam(points_matrix);
        Matrix * ddg_matrix = run_ddg(wam_matrix);
        convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
        Matrix * lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
        print_full_matrix(lnorm_matrix);

        free_matrix_memory(wam_matrix);
        free_matrix_memory(ddg_matrix);
        free_matrix_memory(lnorm_matrix);
    }
    if (strcmp(goal, "spk") == 0) {
        run = 1;
        run_spk(points_matrix, k);
    }

    if (run == 0) {
        printf("Invalid input");
    }

    /***clean all***/
    free_matrix_memory(points_matrix);
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


