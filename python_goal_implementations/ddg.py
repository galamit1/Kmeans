import numpy as np
from .wam import wam_run


def ddg_run(points):
    adjacency_matrix = wam_run(points)
    return ddg_from_adjacency_matrix(adjacency_matrix, len(points))


def ddg_from_adjacency_matrix(adjacency_matrix, n):
    ddg_matrix = np.zeros(adjacency_matrix.shape)
    for i in range(n):
        ddg_matrix[i][i] = 1 / np.sqrt(np.sum(adjacency_matrix[i]))
    return ddg_matrix
