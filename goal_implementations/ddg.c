//
// Created by galam on 26/07/2021.
//

#include "ddg.h"

Matrix * run_dgg(Matrix * wam) {
    Matrix * dgg_matrix = getZerosMatrixSizeN(wam->rows);
    double row_sum;
    for (int i=0; i < wam->rows; i++) {
        row_sum = 0;
        for (int z = i + 1; z < wam->rows; z++) {
            row_sum += wam->cells[i][z];
        }
        dgg_matrix->cells[i][i] = 1 / sqrt(row_sum);
    }
    return dgg_matrix;
}