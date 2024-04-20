from math import sqrt
import os
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/afrq/

To model the Hardy-Weinberg principle, assume that we have a population of N diploid individuals. If an allele is in genetic equilibrium, then because mating is random, we may view the 2N chromosomes as receiving their alleles uniformly. In other words, if there are m dominant alleles, then the probability of a selected chromosome exhibiting the dominant allele is simply p=m2N.

Because the first assumption of genetic equilibrium states that the population is so large as to be ignored, we will assume that N is infinite, so that we only need to concern ourselves with the value of p.

    Given: An array A for which A[k] represents the proportion of homozygous recessive individuals for the k-th Mendelian factor in a diploid population. Assume that the population is in genetic equilibrium for all factors.

    Return: An array B having the same length as A in which B[k] represents the probability that a randomly selected individual carries at least one copy of the recessive allele for the k-th factor.
"""


def verify(result: list[float], solution: list[float]) -> bool:
    correct = len(result) == len(solution) and all(abs(v1-v2) <= 0.001 for v1, v2 in zip(result, solution))
    return correct


def solve(hri: list[float]) -> list[float]:
    # hri = homozygous recessive individuals for the k-th Mendelian factor in a diploid population.
    Bk: list[float] = []
    for ak in Ak:
        p = sqrt(ak)
        q = 1 - p
        prob = ak + 2 * p * q
        Bk.append(prob)
    return Bk


def load_results(path: str) -> list[float]:
    lines = readTextFile(path)
    probs = [float(v) for v in lines[0].split()]

    return probs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    Ak = [float(val) for val in lines[0].split(' ')]

    Bk = solve(Ak)

    solution = load_results(solution_path)
    correct = verify(Bk, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_afrq_1.txt'

    lines = readTextFile(path)
    Ak = [float(val) for val in lines[0].split(' ')]

    Bk = solve(Ak)

    out = " ".join(f"{bk:.4f}" for bk in Bk)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    # correct = solve_and_check(path)
    # print(correct)
