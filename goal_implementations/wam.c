//
// Created by galam on 26/07/2021.
//

#include "wam.h"


Matrix * run_wam(double **points, int num_points, int num_coordinates) {
    Matrix * wam_matrix = getZerosMatrixSizeN(num_points);
    double weight;
    for (int i=0; i < num_points; i++) {
        for (int j = i + 1; j <num_points; j++) {
            weight = calculate_weight(points[i], points[j], num_coordinates);
            wam_matrix->cells[i][j] = weight;
            wam_matrix->cells[j][i] = weight;
        }
    }
    return wam_matrix;
}

double calculate_weight(double * point1, double * point2, int num_coordinates) {
    int weight = 0;
    for (int i = 0; i < num_coordinates; i++) {
        weight += pow(point1[i] - point2[i], 2);
    }
    return exp(-1 * (sqrt(weight) / 2));
}
