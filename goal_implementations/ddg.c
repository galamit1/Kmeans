//
// Created by galam on 26/07/2021.
//

#include "ddg.h"

Matrix * run_ddg(Matrix * wam) {
    Matrix * ddg_matrix = getZerosMatrixSizeN(wam->rows);
    double row_sum;
    for (int i=0; i < wam->rows; i++) {
        row_sum = 0;
        for (int z = 0; z < wam->cols; z++) {
            row_sum += wam->cells[i][z];
        }
        ddg_matrix->cells[i][i] = row_sum;
    }
    return ddg_matrix;
}

void convert_ddg_with_the_pow_of_minus_half(Matrix * ddg_matrix) {
    for (int i=0; i < ddg_matrix->rows; i++) {
        ddg_matrix->cells[i][i] = pow(ddg_matrix->cells[i][i], -0.5);
    }
}