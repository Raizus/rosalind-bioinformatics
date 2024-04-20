
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import math
from BioInfoToolkit.IO.IO import readTextFile


"""
https://rosalind.info/problems/root/

As in the case of unrooted trees, say that we have a fixed collection of n taxa labeling the leaves of a rooted binary tree T. You may like to verify that (by extension of “Counting Phylogenetic Ancestors”) such a tree will contain n-1 internal nodes and 2n-2 total edges. Any edge will still encode a split of taxa; however, the two splits corresponding to the edges incident to the root of T will be equal. We still consider two trees to be equivalent if they have the same splits (which requires that they must also share the same duplicated split to be equal).

Let B(n) represent the total number of distinct rooted binary trees on n labeled taxa.

    Given: A positive integer n (n≤1000).

    Return: The value of B(n) modulo 1,000,000.
"""

def count_rooted_labeled_binary_trees(nleaves: int, mod: int = 1000000) -> int:
    # see https://en.wikipedia.org/wiki/Phylogenetic_tree
    if nleaves >= 2:
        # (2n-3)! / ( (n-2)! * 2^(n-2) )
        count = math.prod(range(2*nleaves-3, n-2, -1))
        count = (count // (2**(nleaves-2))) % mod
        return count

    return 1


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, mod: int) -> int:
    count = count_rooted_labeled_binary_trees(n, mod)
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
    path = f'{cwd}/rosalind_root_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    mod = 10**6

    count = solve(n, mod)

    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    # correct = solve_and_check(path)
    # print(correct)
