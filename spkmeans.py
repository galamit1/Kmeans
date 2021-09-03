import sys

import numpy as np

import mykmeanssp
from python_goal_implementations.wam import wam_run
from python_goal_implementations.ddg import ddg_run
from python_goal_implementations.lnorm import lnorm_run

POINTS_SEPARATOR = '\n'
COORDINATES_SEPARATOR = ','


def get_arguments():
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    with open(file_name, 'r') as f:
        data = f.read()
    points = [[float(num) for num in row.split(COORDINATES_SEPARATOR)] for row in data.split(POINTS_SEPARATOR)]
    if k >= len(points):
        print("Invalid Input!")
        raise Exception
    return k, goal, points

def run_spk(points, k):
    points = np.array(points)
    np.random.seed(0)
    points_num = points.shape[0]
    centroids_indexes = np.random.choice(points_num, 1)
    centroids = points[[centroids_indexes[0]]]
    distances = np.ones(points_num) * float('inf')
    for z in range(1, k):
        for i in range(points_num):
            distances[i] = min(distances[i], np.sum(np.power(np.subtract(points[i], centroids[-1]), 2)))
        probs = np.divide(distances, distances.sum())
        centroids_indexes = np.append(centroids_indexes, np.random.choice(points_num, 1, p=probs), axis=0)
        centroids = np.append(centroids, points[[centroids_indexes[-1]]], axis=0)
    centroids_output = np.array(mykmeanssp.fit(points.tolist(), centroids.tolist(), k, 300, len(points), len(points[0])))
    centroids_output = np.round(centroids_output, decimals=4)
    print(POINTS_SEPARATOR.join([COORDINATES_SEPARATOR.join([str(c) for c in centroid]) for centroid in centroids_output.tolist()]))


def main():
    try:
        k, goal, points = get_arguments()
    except:
        print("Invalid Input!")
        return

    result = []
    if goal == "spk":
        run_spk(points, k)
    elif goal == "wam":
        result = wam_run(points).tolist()
    elif goal == "ddg":
        result = ddg_run(points).tolist()
    elif goal == "lnorm":
        result = lnorm_run(points).tolist()  # we get identity matrix because D*W == zeros !!!
    else:
        print("Invalid goal")
        return

    output = POINTS_SEPARATOR.join(
        [COORDINATES_SEPARATOR.join(["{:.4f}".format(round(i, 4)) for i in c]) for c in result])
    print(output)
    return output


if __name__ == '__main__':
    main()
