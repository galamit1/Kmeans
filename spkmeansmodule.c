#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

/*** Function Declaration ***/
static PyObject* c_kmeans(PyObject* self, PyObject* args);
static PyObject* c_spk(PyObject *self, PyObject *args);
int python_list_of_lists_to_2D_array (PyObject *python_list_of_lists, double **target_array);
Cluster** python_init_k_clusters (int k);
PyMODINIT_FUNC PyInit_myspkmeans(void);


/*******************************/
/*** run functions from goal ***/
/*******************************/

static PyObject* c_spkmeans(PyObject *self, PyObject *args) {
    /*Define variables to receive from user*/
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
        Py_RETURN_NONE; /*In the CPython API, Null is never a valid value for a PyObject* - so it signals an error*/
    }

    /*Verify that data_points & initial_indexes are python lists*/
    if (!PyList_Check(data_points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }

    /*Load points from PyObject into C format of 2D-array points*/
    points = init_points(num_points, num_coordinates);

    if (python_list_of_lists_to_2D_array(data_points, points) != 0) {
        Py_RETURN_NONE;
    }
    if ((points_matrix = malloc(sizeof(Matrix))) == NULL) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }

    points_matrix->rows = num_points;
    points_matrix->cols = num_coordinates;
    points_matrix->cells = points;
    run_functions_according_to_goal(goal, points_matrix, 1);
    Py_RETURN_NONE;
}

/*******************************/
/*** spk function ***/
/*******************************/

static PyObject* c_spk(PyObject *self, PyObject *args) {
    /*Define variables to receive from user*/
    int k;
    int max_iter = 300;
    int num_points;
    int num_coordinates;
    PyObject *data_points;
    PyObject *initial_centroids;
    double **points;
    double **initial_points_for_centroids;
    double *point;
    Cluster **clusters;
    int did_update;
    int i;
    int j;
    PyObject *centroids_output_list;
    Py_ssize_t output_list_len;
    PyObject *single_centroid_list;
    Py_ssize_t output_num_coordinates;
    PyObject *output_coordinate_item;


    /* Parse arguments from Python */
    if((!PyArg_ParseTuple(args, "OOiii", &data_points, &initial_centroids, &k, &num_points, &num_coordinates))) {
        return NULL; /*In the CPython API, Null is never a valid value for a PyObject* - so it signals an error*/
    }


    /*Verify that data_points & initial_indexes are python lists*/
    if (!PyList_Check(data_points) || !PyList_Check(initial_centroids)) {
        return NULL;
    }

    /*Load points from PyObject into C format of 2D-array points*/
    points = init_points(num_points, num_coordinates);

    if (python_list_of_lists_to_2D_array(data_points, points) != 0) {
        return NULL;
    }

    /*Load initial centroids from PyObject into C format of 2D-array -> into clusters  */
    clusters = python_init_k_clusters(k);
    initial_points_for_centroids = init_points(k, num_coordinates);

    if (python_list_of_lists_to_2D_array(initial_centroids, initial_points_for_centroids) != 0) {
        return NULL;
    }

    for (i=0; i<k; i++) {
        clusters[i] = make_cluster(initial_points_for_centroids[i], num_coordinates, i);
    }

    free_points_memory(initial_points_for_centroids, k); /*we don't need this 2D array anymore since initial centroids are now stored in clusters*/

    /***Execute k-means algorithm***/
    did_update = 0;

    for (i=0; i<max_iter; i++) {
        for (j=0; j<num_points; j++) {
            point = points[j];
            find_minimum_distance(clusters, point, k, num_coordinates);
        }
        did_update = update_all_clusters(clusters, k, num_coordinates);
        if (did_update == 1) { /*No update occured - we can break*/
            break;
        }
    }

    free_points_memory(points, num_points);

    /***Convert results to Python Object***/
    output_list_len = k;
    output_num_coordinates = num_coordinates;
    centroids_output_list = PyList_New(output_list_len); /*Create final centroids list*/
    if (centroids_output_list == NULL) {
        return NULL;
    }
    for (i=0; i<k; i++) {
        single_centroid_list = PyList_New(output_num_coordinates); /*Create single centroid list*/
        if (single_centroid_list == NULL) {
            return NULL;
        }
        for (j=0; j<num_coordinates; j++) {
            output_coordinate_item = PyFloat_FromDouble(clusters[i]->curr_centroids[j]);
            if (output_coordinate_item == NULL) {
                return NULL;
            }
            PyList_SET_ITEM(single_centroid_list, j, output_coordinate_item); /*user macro version of PyList_setItem() since there's no previous content*/
        }
        PyList_SET_ITEM(centroids_output_list, i, single_centroid_list);
    }

    free_clusters_memory(clusters, k);

    return centroids_output_list;
}

static PyObject* c_get_t_matrix(PyObject *self, PyObject *args) {
    /*Define variables to receive from user*/
    int k;
    int num_points;
    int num_coordinates;
    PyObject *data_points;
    double **points;
    double *point;
    int i;
    int j;
    Py_ssize_t output_list_len;
    Py_ssize_t output_num_coordinates;
    PyObject *output_list;
    PyObject *single_line;
    PyObject *output_coordinate_item;
    Matrix* points_matrix;
    Matrix* wam_matrix;
    Matrix* ddg_matrix;
    Matrix* lnorm_matrix;
    Matrix* U_matrix;
    Matrix* T_matrix;


    /* Parse arguments from Python */
    if((!PyArg_ParseTuple(args, "Oiii", &data_points, &k, &num_points, &num_coordinates))) {
        return NULL; /*In the CPython API, Null is never a valid value for a PyObject* - so it signals an error*/
    }

    /*Verify that data_points & initial_indexes are python lists*/
    if (!PyList_Check(data_points)) {
        return NULL;
    }

    /*Load points from PyObject into C format of 2D-array points*/
    points = init_points(num_points, num_coordinates);

    if (python_list_of_lists_to_2D_array(data_points, points) != 0) {
        return NULL;
    }

    /* create point matrix */
    points_matrix = get_matrix_from_2D_array(points, num_points, num_coordinates);

    /* run flow to get the T matrix */
    wam_matrix = run_wam(points_matrix);
    ddg_matrix = run_ddg(wam_matrix);
    convert_ddg_with_the_pow_of_minus_half(ddg_matrix);
    lnorm_matrix = run_lnorm(wam_matrix, ddg_matrix);
    U_matrix = run_jacobi(lnorm_matrix, "spk");
    T_matrix = normalize_matrix(U_matrix);

    free_matrix_memory(wam_matrix);
    free_matrix_memory(ddg_matrix);
    free_matrix_memory(lnorm_matrix);
    free_matrix_memory(U_matrix);


    /***Convert results to Python Object***/
    output_list_len = T_matrix->rows;
    output_num_coordinates = T_matrix->cols;
    output_list = PyList_New(output_list_len); /*Create final points list as python object*/
    if (output_list == NULL) {
        return NULL;
    }
    for (i=0; i<T_matrix->rows; i++) {
        single_line = PyList_New(output_num_coordinates); /*Create single line*/
        if (single_line == NULL) {
            return NULL;
        }
        for (j=0; j<T_matrix->cols; j++) {
            output_coordinate_item = PyFloat_FromDouble(T_matrix->cells[i][j]);
            if (output_coordinate_item == NULL) {
                return NULL;
            }
            PyList_SET_ITEM(single_line, j, output_coordinate_item); /*user macro version of PyList_setItem() since there's no previous content*/
        }
        PyList_SET_ITEM(output_list, i, single_line);
    }

    free_matrix_memory(T_matrix);

    return output_list;
}

/*******************************/
/*** Modules setup ***/
/*******************************/

static PyMethodDef _methods[] = {
        {"run_module",
                (PyCFunction) c_spkmeans,
                METH_VARARGS,
                        PyDoc_STR("Run functions according to goal."),
        },
        {"run_spk_module",
                (PyCFunction) c_spk,
                METH_VARARGS,
                PyDoc_STR("Kmeans execution with given initial centroids."),
        },
        {"get_points_for_spk",
            (PyCFunction) c_get_t_matrix,
            METH_VARARGS,
            PyDoc_STR("return the calculated T matrix."),
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

PyMODINIT_FUNC PyInit_myspkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

/*******************************/
/*** Function Implementation ***/
/*******************************/

int python_list_of_lists_to_2D_array (PyObject *python_list_of_lists, double **target_array) {
    Py_ssize_t list_size;
    Py_ssize_t entry_size;
    PyObject *point_item;
    PyObject *coordinate_item;
    int i;
    int j;

    list_size = PyList_Size(python_list_of_lists); /*equivilant to len(_list) in Python*/
    for (i=0; i<list_size; i++) { /*iterate over points*/
        point_item = PyList_GetItem(python_list_of_lists, i);
        if (!PyList_Check(point_item)) {
            return 1;
        }
        entry_size = PyList_Size(point_item);
        for (j=0; j<entry_size; j++) { /*iterate over coordinates of single point*/
            coordinate_item = PyList_GetItem(point_item, j);
            if (!PyFloat_Check(coordinate_item)) {
                return 1;
            }
            target_array[i][j] = PyFloat_AsDouble(coordinate_item);
        }
    }
    return 0;
}

Cluster** python_init_k_clusters (int k) {
    /*
    Recieves pointer to 2D array of points, k, number of coordinates and 2D array of initial indexes,
    Returns new 2D array of Clusters with sufficient memory, initialized with the first k points.
    */
    Cluster **clusters;

    clusters = malloc(sizeof(Cluster*) * k);
    assert (clusters != NULL);

    return clusters;
}