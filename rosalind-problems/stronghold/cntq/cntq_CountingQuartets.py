from BioInfoToolkit.Phylogeny.Phylo import count_quartets
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/cntq/

A quartet AB|CD is consistent with a binary tree T if the quartet can be inferred from one of the splits of T (see “Quartets” for a description of inferring quartets from splits).

Let q(T) denote the total number of quartets that are consistent with T.

    Given: A positive integer n (4≤n≤5000), followed by an unrooted binary tree T in Newick format on n taxa.

    Return: The value of q(T) modulo 1,000,000.
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, mod: int) -> int:
    count = count_quartets(n)
    count = count % mod
    return count


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    nleaves = int(lines[0])
    # newick = lines[1]
    mod = 1000000

    count = solve(nleaves, mod)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(count, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_cntq_1.txt'

    lines = readTextFile(path)
    nleaves = int(lines[0])
    newick = lines[1]
    mod = 1000000

    count = solve(nleaves, mod)
    
    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    correct = solve_and_check(path)
    print(correct)
