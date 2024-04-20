from BioInfoToolkit.Sequences.StringUtils import find_most_frequent_kmers
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1b/

We say that Pattern is a most frequent k-mer in Text if it maximizes Count(Text, Pattern) among all k-mers. For example, "ACTAT" is a most frequent 5-mer in "ACAACTATGCATCACTATCGGGAACTATCCT", and "ATA" is a most frequent 3-mer of "CGATATATCCATAG".

Frequent Words Problem

Find the most frequent k-mers in a string.

    Given: A DNA string Text and an integer k.

    Return: All most frequent k-mers in Text (in any order).
"""

OutputT = list[str]

def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(string: str, k: int) -> OutputT:
    kmers = sorted(find_most_frequent_kmers(string, k))
    return kmers


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    kmers = lines[0].split()
    return kmers


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    string = lines[0]
    k = int(lines[1])
    result = solve(string, k)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1b_1.txt'

    lines = readTextFile(path)
    string = lines[0]
    k = int(lines[1])

    kmers = solve(string, k)

    out = ' '.join(kmers)

    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(out), 'w')

    correct = solve_and_check(path)
    print(correct)

