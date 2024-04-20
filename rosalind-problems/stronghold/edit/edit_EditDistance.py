import BioInfoToolkit.Alignment as Alignment
import os
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/edit/

Given two strings s and t (of possibly different lengths), the edit distance dE(s,t) is the minimum number of edit operations needed to transform s into t, where an edit operation is defined as the substitution, insertion, or deletion of a single symbol.

The latter two operations incorporate the case in which a contiguous interval is inserted into or deleted from a string; such an interval is called a gap. For the purposes of this problem, the insertion or deletion of a gap of length k still counts as k distinct edit operations.

    Given: Two protein strings s and t in FASTA format (each of length at most 1000 aa).

    Return: The edit distance dE(s,t).
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> int:
    dist = Alignment.editDistance(seq1, seq2)
    return dist


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    fasta_dict = read_FASTA(input_path)
    seqs = list(fasta_dict.values())
    seq1 = seqs[0]
    seq2 = seqs[1]

    dist = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(dist, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_edit_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1 = seqs[0]
    seq2 = seqs[1]

    dist = solve(seq1, seq2)

    print(dist)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(dist), 'w')

    # correct = solve_and_check(path)
    # print(correct)
