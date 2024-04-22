from collections import Counter
from BioInfoToolkit.Sets import setFromDiffMultiset, allPairwiseDifferences

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4m/

If A = (a1 = 0, a2, ... , an) is a set of n points on a line segment in increasing order (a1 < a2 < · · · < an), then Δ A denotes the collection of all pairwise differences between points in A. For example, if A = (0, 2, 4, 7), then

    ΔA = (7, 5, 4, 3, 2, 2, 0, 0, 0, 0, 2, 2, 3, 4, 5, 7).

The following problem asks us to reconstruct A from ΔA.

Turnpike Problem
    Given all pairwise distances between points on a line segment, reconstruct the positions of those points.

        Given: A collection of integers L.

        Return: A set A such that ∆A = L.
"""

# https://www.cs.ucf.edu/courses/cap5510/fall2009/res.map/Restrict.Mapping.pdf


def verify(result: list[int], solution: list[int], ms: list[int]) -> bool:
    delta_res = allPairwiseDifferences(result)
    correct = Counter(ms) == Counter(delta_res)
    return correct


def solve(deltaA: list[int]) -> list[int]:
    diffMultiset = [x for x in deltaA if x > 0]
    possibleSets = setFromDiffMultiset(diffMultiset)
    sol1 = sorted(list(possibleSets[0]))
    return sol1


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    set = [int(v) for v in lines[0].split()]
    return set


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    deltaA = [int(a) for a in lines[0].split()]

    possible_set = solve(deltaA)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(possible_set, solution, deltaA)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4m_1.txt'

    lines = readTextFile(path)
    deltaA = [int(a) for a in lines[0].split()]

    result = solve(deltaA)

    out = ' '.join(str(val) for val in result)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

