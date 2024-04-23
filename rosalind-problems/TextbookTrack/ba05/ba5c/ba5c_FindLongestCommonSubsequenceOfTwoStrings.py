from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str, longestCommonSubsequence
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba5c/

Longest Common Subsequence Problem

    Given: Two strings.

    Return: A longest common subsequence of these strings.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    correct = len(result) == len(solution) and is_subsequence_str(
        seq1, result) and is_subsequence_str(seq2, result)
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    string, _ = longestCommonSubsequence(seq1, seq2)
    return string


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lcs = lines[0]
    return lcs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    seq1 = lines[0]
    seq2 = lines[1]

    result = solve(seq1, seq2)

    solution = load_results(solution_path)

    correct = verify(result, solution, seq1, seq2)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba5c_1.txt'

    lines = readTextFile(path)
    seq1 = lines[0]
    seq2 = lines[1]

    string = solve(seq1, seq2)

    out = string
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

