import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
A partial permutation is an ordering of only k objects taken from a collection containing n objects (i.e., k≤n). For example, one partial permutation of three of the first eight positive integers is given by (5,7,2).

The statistic P(n,k) counts the total number of partial permutations of k objects that can be formed from a collection of n objects. Note that P(n,n) is just the number of permutations of n objects, which we found to be equal to n!=n(n-1)(n-2)⋯(3)(2) in “Enumerating Gene Orders”.

    Given: Positive integers n and k such that 100≥n>0 and 10≥k>0.

    Return: The total number of partial permutations P(n,k), modulo 1,000,000.
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, k: int, mod: int) -> int:
    Pn_k = 1
    for i in range(n, n-k, -1):
        Pn_k = (Pn_k * i) % mod
    return Pn_k


def load_results(path: str) -> int:
    lines = readTextFile(path)
    perms = int(lines[0])
    return perms


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    mod = 1000000
    n, k = [int(v) for v in lines[0].split()]
    result = solve(n, k, mod)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_pper_1.txt'

    lines = readTextFile(path)
    mod = 1000000
    n, k = [int(v) for v in lines[0].split()]
    Pn_k = solve(n,k,mod)

    print(Pn_k)

    result_path = result_path_from_input_path(path)

    writeTextFile(result_path, str(Pn_k), 'w')

    correct = solve_and_check(path)
    print(correct)
