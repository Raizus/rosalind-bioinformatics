import math
from scipy.stats import binom
import numpy as np
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/foun/

    Given: Two positive integers N and m, followed by an array A containing k integers between 0 and 2N. A[j] represents the number of recessive alleles for the j-th factor in a population of N diploid individuals.

    Return: An m x k matrix B for which Bi,j represents the common logarithm of the probability that after i generations, no copies of the recessive allele for the j-th factor will remain in the population. Apply the Wright-Fisher model.
"""


def wrightFisherNextGen(p: list[float], n: int):
    assert len(p) == n+1
    return [sum([float(binom.pmf(r, n, i/n))*p[i] for i in range(n+1)])
            for r in range(1+n)]


def wrightFisherModel(n: int, m: int, g: int) -> list[float]:
    """Cumputes the probability mass function pmf(x), where P(X=k), k in [0,m], is the probability that there are k alleles in the population after g generations

    Args:
        n (int): Number of chromosomes in the population
        m (int): Number of alleles in the population
        g (int): Number of gerenations

    Returns:
        list[float]: _description_
    """
    ms: list[float] = [1.0 if i == m else 0 for i in range(n+1)]
    for _ in range(g):
        ms = wrightFisherNextGen(ms, n)
    return ms


def founderMatrix(n: int, g: int, Ak: list[int]):
    num_factors = len(Ak)
    mat = np.zeros((g, num_factors))
    for j, m_recessive in enumerate(Ak):
        ms: list[float] = [1.0 if i == m_recessive else 0 for i in range(n+1)]
        for i in range(g):
            ms = wrightFisherNextGen(ms, n)
            mat[i, j] = math.log10(ms[0])
    return mat


OutputT = list[list[float]]

def verify(result: OutputT, solution: OutputT) -> bool:
    if len(result) != len(solution):
        return False
    for line1, line2 in zip(result, solution):
        if len(line1) != len(line2) or not all(abs(v1-v2) <= 0.001 for v1, v2 in zip(line1, line2)):
            return False
    return True


def solve(N: int, m: int, Ak: list[int]) -> OutputT:
    mat: list[list[float]] = founderMatrix(2*N, m, Ak).tolist()
    return mat


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    mat: OutputT = [[float(v) for v in line.split()] for line in lines]
    return mat


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    N, m = [int(v) for v in lines[0].split()]
    Ak = [int(v) for v in lines[1].split()]

    result = solve(N, m, Ak)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_foun_1.txt'

    lines = readTextFile(path)
    N, m = [int(v) for v in lines[0].split()]
    Ak = [int(v) for v in lines[1].split()]

    mat = founderMatrix(2*N, m, Ak)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for row in mat:
        out = " ".join(f"{v:.8f}" for v in row)
        print(out)
        writeTextFile(result_path, out, 'a')


    correct = solve_and_check(path)
    print(correct)

