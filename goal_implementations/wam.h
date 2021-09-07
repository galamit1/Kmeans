//
// Created by galam on 26/07/2021.
//
#ifndef SPK_WAM_H
#define SPK_WAM_H

#include "../utils.h"

Matrix * run_wam(Matrix * points);
double calculate_weight(double * point1, double * point2, int num_coordinates);

#endif //SPK_WAM_H