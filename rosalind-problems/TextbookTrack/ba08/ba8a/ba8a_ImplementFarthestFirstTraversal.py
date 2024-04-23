from BioInfoToolkit.Clustering import farthestFirstTravel
import numpy as np
import numpy.typing as npt
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba8a/

To think about clustering as dividing data points Data into k clusters, we will try to select a set Centers of k points that will serve as the centers of these clusters. We would like to choose Centers so that they minimize some distance function between Centers and Data over all possible choices of centers. But how should this distance function be defined?

First, we define the Euclidean distance between points v = (v1, ... , vm) and w = (w1, ... , wm) in m-dimensional space, denoted d(v, w), as the length of the line segment connecting these points,

    d(v,w)=√ [∑_(i=1)^m (vi-wi)^2]

Next, given a point DataPoint in multi-dimensional space and a set of k points Centers, we define the distance from DataPoint to Centers, denoted d(DataPoint, Centers), as the Euclidean distance from DataPoint to its closest center,

    d(DataPoint,Centers) = minall points x from Centersd(DataPoint, x).

The length of the segments in the figure below correspond to d(DataPoint, Centers) for each point DataPoint.

We now define the distance between all data points Data and centers Centers. This distance, denoted MaxDistance(Data, Centers), is the maximum of d(DataPoint, Centers) among all data points DataPoint,

    MaxDistance(Data, Centers) = maxall points DataPoint from Data d(DataPoints,Centers).

We can now formulate a well-defined clustering problem.

k-Center Clustering Problem: Given a set of data points, find k centers minimizing the maximum distance between these data points and centers.
    Input: A set of points Data and an integer k.
    Output: A set Centers of k centers that minimize the distance MaxDistance(DataPoints, Centers) over all possible choices of k centers.

Although the k-Center Clustering Problem is easy to state, it is NP-Hard. The Farthest First Traversal heuristic, whose pseudocode is shown below, selects centers from the points in Data (instead of from all possible points in m-dimensional space). It begins by selecting an arbitrary point in Data as the first center and iteratively adds a new center as the point in Data that is farthest from the centers chosen so far, with ties broken arbitrarily.

    FarthestFirstTraversal(Data, k) 
        Centers ← the set consisting of a single randomly chosen point from Data
            while |Centers| < k
                DataPoint ← the point in Data maximizing d(DataPoint, Centers) 
                add DataPoint to Centers 
        return Centers 

Implement FarthestFirstTraversal

    Given: Integers k and m followed by a set of points Data in m-dimensional space.

    Return: A set Centers consisting of k points (centers) resulting from applying FarthestFirstTraversal(Data, k), where the first point from Data is chosen as the first center to initialize the algorithm.
"""

OutputT = list[list[float]]


def verify(result: OutputT, solution: OutputT) -> bool:
    def point_matches(p1: list[float], p2: list[float]):
       return len(p1) == len(p2) and all(abs(a-b) <= 0.001 for a, b in zip(p1, p2))

    if len(result) != len(solution):
        return False
    if not point_matches(result[0], solution[0]):
        return False

    # TODO: check the rest

    return True


def solve(data: npt.NDArray, k: int) -> OutputT:
    centers = farthestFirstTravel(data, k)
    centers = centers.tolist()
    return centers


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    centers: OutputT = [[float(v) for v in line.split()] for line in lines]
    return centers


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, m = [int(v) for v in lines[0].split()]
    data = [[float(v) for v in line.split()] for line in lines[1:]]
    data = np.array(data)

    result = solve(data, k)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba8a_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, m = [int(v) for v in lines[0].split()]
    data = [[float(v) for v in line.split()] for line in lines[1:]]
    data = np.array(data)

    centers = solve(data, k)


    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for center in centers:
        out = ' '.join(str(v) for v in center)
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
