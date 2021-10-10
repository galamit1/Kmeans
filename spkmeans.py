import sys

import numpy as np

import myspkmeans

POINTS_SEPARATOR = '\n'
COORDINATES_SEPARATOR = ','
normalize_num = lambda num: f"{num * -1:.4f}" if -0.00005 < num < 0 else f"{num:.4f}"


def get_arguments():
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    with open(file_name, 'r') as f:
        data = f.read()
    rows = [row for row in data.split(POINTS_SEPARATOR) if row != ""]
    points = [[float(num.strip()) for num in row.split(COORDINATES_SEPARATOR)] for row in rows]
    return k, goal, points

def run_spk(points, k):
    if k >= len(points):
        print("Invalid Input!")
        raise Exception
    points = np.array(points)
    np.random.seed(0)
    points_num = points.shape[0]
    centroids_indexes = np.random.choice(points_num, 1)
    centroids = points[[centroids_indexes[0]]]
    distances = np.ones(points_num) * float('inf')
    for _ in range(1, k):
        for i in range(points_num):
            distances[i] = min(distances[i], np.sum(np.power(np.subtract(points[i], centroids[-1]), 2)))
        probs = np.divide(distances, distances.sum())
        centroids_indexes = np.append(centroids_indexes, np.random.choice(points_num, 1, p=probs), axis=0)
        centroids = np.append(centroids, points[[centroids_indexes[-1]]], axis=0)
    print(COORDINATES_SEPARATOR.join([str(c) for c in centroids_indexes]))
    centroids_output = myspkmeans.run_spk_module(points.tolist(), centroids.tolist(), k, len(points), len(points[0]))
    print(POINTS_SEPARATOR.join(
        [COORDINATES_SEPARATOR.join([normalize_num(c) for c in centroid]) for centroid in centroids_output]))


def main():
    k, goal, points = get_arguments()
    if goal == "spk":
        t_matrix = myspkmeans.get_points_for_spk(points, k, len(points), len(points[0]))
        run_spk(t_matrix, len(t_matrix[0]))
    else:
        myspkmeans.run_module(points, goal, len(points), len(points[0]))


if __name__ == '__main__':
    main()
