from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str, shortestCommonSupersequence
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/scsp/

A string s is a supersequence of another string t if s contains t as a subsequence.

A common supersequence of strings s and t is a string that serves as a supersequence of both s and t. For example, "GACCTAGGAACTC" serves as a common supersequence of "ACGTC" and "ATAT". A shortest common supersequence of s and t is a supersequence for which there does not exist a shorter common supersequence. Continuing our example, "ACGTACT" is a shortest common supersequence of "ACGTC" and "ATAT".

Given: Two DNA strings s and t.

Return: A shortest common supersequence of s and t. If multiple solutions exist, you may output any one.
"""


def verify(result: str, solution: str, seq1: str, seq2: str) -> bool:
    correct = len(result) == len(solution) and is_subsequence_str(
        result, seq1) and is_subsequence_str(result, seq2)
    return correct


def solve(s: str, t: str) -> str:
    scs = shortestCommonSupersequence(s, t)
    return scs


def load_results(path: str) -> str:
    lines = readTextFile(path)
    scs = lines[0]
    return scs


def solve_and_check(input_path: str) -> bool:
    seqs = readTextFile(input_path)
    s, t = seqs[0], seqs[1]

    scs = solve(s,t)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(scs, solution, s, t)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_scsp_1.txt'

    seqs = readTextFile(path)
    s,t = seqs[0], seqs[1]

    scs = shortestCommonSupersequence(s, t)

    print(scs)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, scs, 'w')

    correct = solve_and_check(path)
    print(correct)

