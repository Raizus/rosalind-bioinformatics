
import math
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/aspc/

In “Counting Subsets”, we saw that the total number of subsets of a set S containing n elements is equal to 2n.

However, if we intend to count the total number of subsets of S having a fixed size k, then we use the combination statistic C(n,k), also written (nk).

    Given: Positive integers n and m with 0≤m≤n≤2000.

    Return: The sum of combinations C(n,k) for all k satisfying m≤k≤n, modulo 1,000,000. In shorthand, ∑nk=m(nk).
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, m: int, mod: int) -> int:
    count = 0
    for k in range(m, n+1):
        count = (count + math.comb(n, k)) % mod
    return count


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n, m = [int(v) for v in lines[0].split()]
    mod = 1000000

    count = solve(n, m, mod)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(count, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_aspc_1.txt'

    lines = readTextFile(path)
    n, m = [int(v) for v in lines[0].split()]
    mod = 1000000

    count = solve(n, m, mod)

    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    # correct = solve_and_check(path)
    # print(correct)
