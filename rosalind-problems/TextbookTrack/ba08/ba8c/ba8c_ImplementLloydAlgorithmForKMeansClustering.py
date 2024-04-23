from BioInfoToolkit.Clustering import LloydAlgorithm
import numpy as np
import numpy.typing as npt
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba8c/

The Lloyd algorithm is one of the most popular clustering heuristics for the k-Means Clustering Problem. It first chooses k arbitrary points Centers from Data as centers and then iteratively performs the following two steps:

    Centers to Clusters: After centers have been selected, assign each data point to the cluster corresponding to its nearest center; ties are broken arbitrarily. 
    Clusters to Centers: After data points have been assigned to clusters, assign each cluster's center of gravity to be the cluster's new center. 

We say that the Lloyd algorithm has converged if the centers (and therefore their clusters) stop changing between iterations. 

Implement the Lloyd algorithm

    Given: Integers k and m followed by a set of points Data in m-dimensional space.

    Return: A set Centers consisting of k points (centers) resulting from applying the Lloyd algorithm to Data and Centers, where the first k points from Data are selected as the first k centers. 
"""

OutputT = list[list[float]]


def verify(result: OutputT, solution: OutputT) -> bool:
    def point_matches(p1: list[float], p2: list[float]):
       return len(p1) == len(p2) and all(abs(a-b) <= 0.001 for a, b in zip(p1, p2))

    idxs: list[int] = []
    for c1 in result:
        idx = [i for i, c2 in enumerate(solution) if point_matches(c1, c2)]
        if len(idx):
            idxs.extend(idx)

    idxs = sorted(idxs)
    correct = all(i == j for i,j in zip(idxs, range(len(solution)))) and len(idxs) == len(solution)

    return correct


def solve(data: npt.NDArray, k: int) -> OutputT:
    centers = LloydAlgorithm(data, k, 'first_k')
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
    path = f'{cwd}/rosalind_ba8c_1.txt'

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
