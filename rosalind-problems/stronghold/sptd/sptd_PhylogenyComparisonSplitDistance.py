from BioInfoToolkit.Phylogeny.Phylo import PhyloTree, split_distance
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/sptd/

Define the split distance between two unrooted binary trees as the number of nontrivial splits contained in one tree but not the other.

Formally, if s(T1,T2) denotes the number of nontrivial splits shared by unrooted binary trees T1 and T2, Then their split distance is dsplit(T1,T2)=2(n-3)-2s(T1,T2).

    Given: A collection of at most 3,000 species taxa and two unrooted binary trees T1 and T2 on these taxa in Newick format.

    Return: The split distance dsplit(T1,T2).
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(taxa: list[str], newick1: str, newick2: str) -> int:
    phylo1 = PhyloTree(newick1)
    phylo2 = PhyloTree(newick2)

    dist = split_distance(taxa, phylo1, phylo2)
    return dist


def load_results(path: str) -> int:
    lines = readTextFile(path)
    dist = int(lines[0])
    return dist


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    taxa = lines[0].split(' ')
    newick1 = lines[1]
    newick2 = lines[2]

    dist = solve(taxa, newick1, newick2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(dist, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_sptd_1.txt'

    lines = readTextFile(path)
    taxa = lines[0].split(' ')
    newick1 = lines[1]
    newick2 = lines[2]

    dist = solve(taxa, newick1, newick2)
    print(dist)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(dist), 'w')

    correct = solve_and_check(path)
    print(correct)

