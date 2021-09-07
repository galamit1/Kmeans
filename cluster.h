//
// Created by galam on 07/09/2021.
//

#ifndef SPK_CLUSTER_H
#define SPK_CLUSTER_H


typedef struct cluster {
    int cluster_size;
    int cluster_index;
    double* prev_centroids;
    double* curr_centroids;
    double* sum;
} Cluster;


#endif //SPK_CLUSTER_H
