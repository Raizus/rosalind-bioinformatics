
import math
from BioInfoToolkit.IO.IO import readTextFile

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/cunr/

Two unrooted binary trees T1 and T2 having the same n labeled leaves are considered to be equivalent if there is some assignment of labels to the internal nodes of T1 and T2 so that the adjacency lists of the two trees coincide. As a result, note that T1 and T2 must have the same splits; conversely, if the two trees do not have the same splits, then they are considered distinct.

Let b(n) denote the total number of distinct unrooted binary trees having n labeled leaves.

    Given: A positive integer n (n≤1000).

    Return: The value of b(n) modulo 1,000,000.
"""


def count_unrooted_labeled_binary_trees(nleaves: int, mod: int = 1000000) -> int:
    # see https://en.wikipedia.org/wiki/Phylogenetic_tree
    if nleaves >= 3:
        # (2n-5)! / ( (n-3)! * 2^(n-3) )
        count = math.prod(range(2*nleaves-5, n-3, -1))
        count = (count // (2**(nleaves-3))) % mod
        return count

    return 1



def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, mod: int) -> int:
    count = count_unrooted_labeled_binary_trees(n, mod)
    return count


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = int(lines[0])
    mod = 10**6

    count = solve(n, mod)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(count, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_cunr_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    mod = 10**6

    count = solve(n, mod)

    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    correct = solve_and_check(path)
    print(correct)
