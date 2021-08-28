#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "utils.h"
#include "goal_implementations/wam.h"
#include "goal_implementations/ddg.h"
#include "goal_implementations/lnorm.h"
#include "goal_implementations/spk.h"

int main(int argc, char **argv);
int main1(int argc, char **argv);
int get_num_coordinates(char* sinlge_line);
double** init_points (int num_points, int num_coordinates);
double* convert_line_to_point (double* point, char* line);
void free_points_memory (double** points, int num_points);
