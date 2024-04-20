from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str, longestCommonSubsequence
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/lcsq/

A string u is a common subsequence of strings s and t if the symbols of u appear in order as a subsequence of both s and t. For example, "ACTG" is a common subsequence of "AACCTTGG" and "ACACTGTGA".

Analogously to the definition of longest common substring, u is a longest common subsequence of s and t if there does not exist a longer common subsequence of the two strings. Continuing our above example, "ACCTTG" is a longest common subsequence of "AACCTTGG" and "ACACTGTGA", as is "AACTGG".

    Given: Two DNA strings s and t (each having length at most 1 kbp) in FASTA format.

    Return: A longest common subsequence of s and t. (If more than one solution exists, you may return any one.)
"""


def verify(result: str, solution: str, seq1: str, seq2: str) -> bool:
    correct = len(result) == len(solution) and is_subsequence_str(
        seq1, result) and is_subsequence_str(seq2, result)
    return correct


def solve(seq1: str, seq2: str) -> str:
    lcs, _ = longestCommonSubsequence(seq1, seq2)
    return lcs


def load_results(path: str) -> str:
    lines = readTextFile(path)
    lcs = lines[0]
    return lcs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1, seq2 = seqs[0], seqs[1]

    result = solve(seq1, seq2)

    solution = load_results(solution_path)

    correct = verify(result, solution, seq1, seq2)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_lcsq_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1, seq2 = seqs[0], seqs[1]

    lcs, _ = longestCommonSubsequence(seq1, seq2)

    print(lcs)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, lcs, 'w')

    # correct = solve_and_check(path)
    # print(correct)

