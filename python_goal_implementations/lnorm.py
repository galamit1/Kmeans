import numpy as np
from .wam import wam_run
from .ddg import ddg_from_adjacency_matrix


def lnorm_run(points):
    adjacency_matrix = wam_run(points)
    ddg_matrix = ddg_from_adjacency_matrix(adjacency_matrix, len(points))
    unit_matrix = np.identity(len(points))
    return unit_matrix - np.multiply(np.multiply(ddg_matrix, adjacency_matrix), ddg_matrix)





