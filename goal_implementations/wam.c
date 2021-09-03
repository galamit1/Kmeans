//
// Created by galam on 26/07/2021.
//

#include "wam.h"


Matrix * run_wam(Matrix * points) {
    Matrix * wam_matrix = getZerosMatrixSizeN(points->rows);
    double weight;
    for (int i=0; i < points->rows; i++) {
        for (int j = i + 1; j < points->rows; j++) {
            weight = calculate_weight(points->cells[i], points->cells[j], points->cols);
            wam_matrix->cells[i][j] = weight;
            wam_matrix->cells[j][i] = weight;
        }
    }
    return wam_matrix;
}

double calculate_weight(double * point1, double * point2, int num_coordinates) {
    double weight = 0;
    for (int i = 0; i < num_coordinates; i++) {
        weight += pow(point1[i] - point2[i], 2);
    }
    return exp(-1 * (sqrt(weight) / 2));
}
