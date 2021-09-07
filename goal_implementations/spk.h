#ifndef SPK_SPK_H
#define SPK_SPK_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define EPSILON 0.001
#define MAX_ITER 500

#include "../cluster.h"
#include "../utils.h"
#include "../kmeans.h"



/*** Function Declaration ***/
void run_spk(Matrix * points, int k);
Cluster** init_k_clusters (Matrix * points, int k);

#endif //SPK_SPK_H