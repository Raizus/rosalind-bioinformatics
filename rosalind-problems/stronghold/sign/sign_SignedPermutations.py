
from itertools import permutations, product
from math import factorial
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
A signed permutation of length n is some ordering of the positive integers {1,2,…,n} in which each integer is then provided with either a positive or negative sign (for the sake of simplicity, we omit the positive sign). For example, π=(5,-3,-2,1,4) is a signed permutation of length 5.

    Given: A positive integer n≤6.

    Return: The total number of signed permutations of length n, followed by a list of all such permutations (you may list the signed permutations in any order).
"""


def verify(result: tuple[int, list[list[int]]], solution: tuple[int, list[list[int]]]) -> bool:
    count_res, perms_res = result
    count_sol, perms_sol = solution
    c1 = count_res == count_sol
    c2 = all(l1 in perms_sol for l1 in perms_res)
    return c1 and c2


def solve(n: int) -> tuple[int, list[list[int]]]:
    count = factorial(n) * n**2
    perms: list[list[int]] = []
    for p1 in permutations(range(1, n+1)):
        for p2 in product([-1, 1], repeat=n):
            perm = [_p1*_p2 for _p1, _p2 in zip(p1, p2)]
            perms.append(perm)

    return count, perms


def load_results(path: str) -> tuple[int, list[list[int]]]:
    lines = readTextFile(path)
    count = int(lines[0])
    perms = [[int(v) for v in line.split()] for line in lines[1:]]
    return count, perms


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    n = int(lines[0])

    count, perms = solve(n)

    solution = load_results(solution_path)

    correct = verify((count, perms), solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_sign_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])

    count, perms = solve(n)

    result_path = result_path_from_input_path(path)

    print(count)
    writeTextFile(result_path, str(count), 'w')
    for perm in perms:
        out = ' '.join(str(p) for p in perm)
        writeTextFile(result_path, out, 'a')
        print(out)

    correct = solve_and_check(path)
    print(correct)
