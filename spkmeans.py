import sys

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


def main():
    try:
        k, goal, points = get_arguments()
    except:
        print("Invalid Input!")
        return

    result = []
    if goal == "spk":
        pass
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
