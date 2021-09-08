//
// Created by galam on 08/09/2021.
//
#include "spkmeans.c"

#ifndef SPK_SPKMEANSMODULE_H
#define SPK_SPKMEANSMODULE_H

#define KMEANS_EPSILON 0.001

/*** Function Declaration ***/
static PyObject* c_kmeans(PyObject* self, PyObject* args);
static PyObject* c_spk(PyObject *self, PyObject *args);
int python_list_of_lists_to_2D_array (PyObject *python_list_of_lists, double **target_array);
Cluster** python_init_k_clusters (int k);
PyMODINIT_FUNC PyInit_myspkmeans(void);


#endif //SPK_SPKMEANSMODULE_H
