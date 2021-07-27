import numpy as np


def wam_run(points):
    points_number = len(points)
    points = np.array(points)
    adjacency_matrix = np.zeros((points_number, points_number))
    for i in range(points_number):
        for j in range(i + 1, points_number):
            distance = np.sqrt(np.sum(np.power(points[i] - points[j], 2)))
            weight = pow(np.e, -0.5 * distance)
            adjacency_matrix[i][j] = weight
            adjacency_matrix[j][i] = weight
    return adjacency_matrix

