from BioInfoToolkit.Clustering import hierarquicalClustering
import numpy as np
import numpy.typing as npt
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba8d/

The soft k-means clustering algorithm starts from randomly chosen centers and iterates the following two steps:

    - Centers to Soft Clusters (E-step): After centers have been selected, assign each data point a “responsibility” value for each cluster, where higher values correspond to stronger cluster membership.
    - Soft Clusters to Centers (M-step): After data points have been assigned to soft clusters, compute new centers.

We begin with the “Centers to Soft Clusters” step. If we think about the centers as stars and the data points as planets, then the closer a point is to a center, the stronger that center’s “pull” should be on the point. Given k centers Centers = (x1, ..., xk) and n points Data = (Data1, ... , Datan), we therefore need to construct a k × n responsibility matrix HiddenMatrix for which HiddenMatrixi,j is the pull of center i on data point j. This pull can be computed according to the Newtonian inverse-square law of gravitation,

    HiddenMatrixi,j=1/d(Dataj,xi)2∑all centers xi1/d(Dataj,xi)2

Unfortunately for Newton fans, the following partition function from statistical physics often works better in practice:

    HiddenMatrixi,j=e-β⋅d(Dataj,xi)∑all centers xie-β⋅d(Dataj,xi)

In this formula, e is the base of the natural logarithm (e ≈ 2.72), and β is a parameter reflecting the amount of flexibility in our soft assignment and called — appropriately enough — the stiffness parameter.

In soft k-means clustering, if we let HiddenMatrixi denote the i-th row of HiddenMatrix, then we can update center xi using an analogue of the above formulas. Specifically, we will define the j-th coordinate of center xi, denoted xi, j, as

    xi,j=HiddenMatrixi⋅DatajHiddenMatrixi⋅1→

Here, Dataj is the n-dimensional vector holding the j-th coordinates of the n points in Data.

The updated center xi is called a weighted center of gravity of the points Data.

Implement the Soft k-Means Clustering Algorithm

    Given: Integers k and m, followed by a stiffness parameter β, followed by a set of points Data in m-dimensional space.

    Return: A set Centers consisting of k points (centers) resulting from applying the soft k-means clustering algorithm. Select the first k points from Data as the first centers for the algorithm and run the algorithm for 100 steps. Results should be accurate up to three decimal places.
"""

OutputT = list[set[int]]


def verify(result: OutputT, solution: OutputT) -> bool:
    idxs: list[int] = []
    for c1 in result:
        idx = [i for i, c2 in enumerate( solution) if c1 == c2]
        if len(idx):
            idxs.extend(idx)

    idxs = sorted(idxs)
    correct = all(i == j for i, j in zip(idxs, range(
        len(solution)))) and len(idxs) == len(solution)

    return correct


def solve(data: npt.NDArray) -> OutputT:
    new_clusters = hierarquicalClustering(data)
    new_clusters: OutputT = [set([v+1 for v in c]) for c in new_clusters]
    return new_clusters


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    clusters: OutputT = []
    for line in lines:
        c = set([int(v) for v in line.split()])
        clusters.append(c)
    return clusters


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    n = int(lines[0])

    data = [[float(v) for v in line.split()] for line in lines[1:]]
    data = np.array(data)

    result = solve(data)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba8e_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    n = int(lines[0])

    data = [[float(v) for v in line.split()] for line in lines[1:]]
    data = np.array(data)

    centers = solve(data)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for center in centers:
        out = ' '.join(str(v) for v in center)
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
