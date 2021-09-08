import sys

import numpy as np

import myspkmeans

POINTS_SEPARATOR = '\n'
COORDINATES_SEPARATOR = ','


def get_arguments():
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    with open(file_name, 'r') as f:
        data = f.read()
    points = [[float(num.strip()) for num in row.split(COORDINATES_SEPARATOR)] for row in data.split(POINTS_SEPARATOR)]
    return k, goal, points

def run_spk(points, k, goal):
    if k >= len(points):
        print("Invalid Input!")
        raise Exception
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
    centroids_output = np.array(myspkmeans.run_spk_module(points.tolist(), centroids.tolist(), k, len(points), len(points[0])))
    centroids_output = np.round(centroids_output, decimals=4)
    print(POINTS_SEPARATOR.join([COORDINATES_SEPARATOR.join([str(c) for c in centroid]) for centroid in centroids_output.tolist()]))


def main():
    k, goal, points = get_arguments()
    if goal == "spk" and k != 0: # I'm not sure about the k != 0
        run_spk(points, k, goal)
    else:
        myspkmeans.run_module(points, goal, len(points), len(points[0]))

if __name__ == '__main__':
    main()
