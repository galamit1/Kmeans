//
// Created by galam on 26/07/2021.
//

#include "spk.h"


void run_spk(Matrix * points, int k) {
    /*argc = number of parameters (including program name), **argv = string array of parameters.*/

    /***Variables declarations***/
    double *point;
    Cluster **clusters;
    int did_update;


    /***Execute k-means algorithm***/
    clusters = init_k_clusters(points, k);
    did_update = 0;

    for (int i=0; i<MAX_ITER; i++) {
        for (int j=0; j<points->rows; j++) {
            point = points->cells[j];
            find_minimum_distance(clusters, point, k, points->cols);
        }
        did_update = update_all_clusters(clusters, k, points->cols);
        if (did_update == 1) { /*No update occured - we can break*/
            break;
        }
    }

    /***Print results***/
    for (int i=0; i<k; i++) {
        for (int j=0; j<points->cols; j++) {
            printf("%.4f", clusters[i]->curr_centroids[j]);
            if (j != (points->cols -1)) {
                printf(",");
            }
        }
        printf("\n");
    }
    free_clusters_memory(clusters, k);
}


Cluster** init_k_clusters (Matrix * points, int k) {
    /*
    Recieves pointer to 2D array of points, k and number of coordinates,
    Returns new 2D array of Clusters with sufficient memory, initialized with the first k points.
    */
    Cluster **clusters;

    clusters = malloc(sizeof(Cluster*) * k);
    assert (clusters != NULL);

    for (int i=0; i<k; i++) {
        clusters[i] = make_cluster(points->cells[i], points->cols, i);
    }

    return clusters;
}
