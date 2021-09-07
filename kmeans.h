//
// Created by galam on 07/09/2021.
//

#ifndef SPK_KMEANS_H
#define SPK_KMEANS_H

#include "cluster.h"
#define EPSILON 0.001

/*** Function Declaration ***/
static PyObject* c_kmeans(PyObject *self, PyObject *args);
double** init_points (int num_points, int num_coordinates);
int python_list_of_lists_to_2D_array (PyObject *python_list_of_lists, double **target_array);
Cluster* make_cluster (const double* point, const int num_coordinates, const int index);
Cluster** python_init_k_clusters (int k);
void add_point_to_cluster (Cluster* cluster, const double* point, int num_coordinates);
void find_minimum_distance (Cluster** clusters, double* point, int k, int num_coordinates);
double update_cluster_centroids_and_sum (Cluster* cluster, int num_coordinates);
int update_all_clusters (Cluster** clusters, int k, int num_coordinates);
void free_clusters_memory (Cluster** clusters, int k);
void free_points_memory (double** points, int num_points);
double get_distance (Cluster* cluster, const double* point, int num_coordinates);

#endif //SPK_KMEANS_H
