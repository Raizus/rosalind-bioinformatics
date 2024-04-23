from BioInfoToolkit.Clustering import squaredErrorDistortion
import numpy as np
import numpy.typing as npt
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba8b/

FarthestFirstTraversal (which we introduced in “Implement FarthestFirstTraversal”) is fast, and its solution approximates the optimal solution of the k-Center Clustering Problem; however, this algorithm is rarely used for gene expression analysis. In k-Center Clustering, we selected Centers so that these points would minimize MaxDistance(Data, Centers), the maximum distance between any point in Data and its nearest center. But biologists are usually interested in analyzing typical rather than maximum deviations, since the latter may correspond to outliers representing experimental errors.

To address limitations of MaxDistance, we will introduce a new scoring function. Given a set Data of n data points and a set Centers of k centers, the squared error distortion of Data and Centers, denoted Distortion(Data, Centers), is defined as the mean squared distance from each data point to its nearest center,

    Distortion(Data,Centers) = (1/n) ∑all points DataPoint in Datad(DataPoint, Centers)2.

Squared Error Distortion Problem

    Given: Integers k and m, followed by a set of centers Centers and a set of points Data.

    Return: The squared error distortion Distortion(Data, Centers).
"""

OutputT = float


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = abs(result-solution) <=0.001
    return correct


def solve(data: npt.NDArray, centers: npt.NDArray) -> OutputT:
    distortion = squaredErrorDistortion(data, centers)
    return float(distortion)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    distortion = float(lines[0])

    return distortion


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, m = [int(v) for v in lines[0].split()]

    centers = [[float(v) for v in line.split()] for line in lines[1:1+k]]
    centers = np.array(centers)

    data = [[float(v) for v in line.split()] for line in lines[1+k+1:]]
    data = np.array(data)

    result = solve(data, centers)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba8b_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, m = [int(v) for v in lines[0].split()]

    centers = [[float(v) for v in line.split()] for line in lines[1:1+k]]
    centers = np.array(centers)

    data = [[float(v) for v in line.split()] for line in lines[1+k+1:]]
    data = np.array(data)

    distortion = solve(data, centers)

    out = str(distortion)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

