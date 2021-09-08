#ifndef SPK_SPKMEANS_H
#define SPK_SPKMEANS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

//#include "utils.h"
#include "kmeans.h"
#include "utils.h"
#include "goal_implementations/wam.h"
#include "goal_implementations/ddg.h"
#include "goal_implementations/lnorm.h"
#include "goal_implementations/spk.h"
#include "goal_implementations/jacobi.h"



int main(int argc, char **argv);
void run_functions_according_to_goal(char * goal, Matrix * points_matrix, int k);


#endif //SPK_SPKMEANS_H