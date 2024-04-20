
from itertools import permutations
from math import factorial
from BioInfoToolkit.IO import writeTextFile
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
A permutation of length n is an ordering of the positive integers {1,2,…,n}. For example, π=(5,3,2,1,4) is a permutation of length 5.

    Given: A positive integer n≤7.

    Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int) -> int:
    perm = factorial(n)
    return perm


def load_results(path: str) -> int:
    lines = readTextFile(path)
    perms = int(lines[0])
    return perms


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    n = int(lines[0])

    result = solve(n)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_perm_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])

    perm = factorial(n)
    values = list(range(1, n+1))
    print(perm)

    result_path = result_path_from_input_path(path)

    writeTextFile(result_path, str(perm), 'w')
    for p in permutations(values, n):
        out = " ".join(str(v) for v in p)
        print(out)
        writeTextFile(result_path, out, 'a')
