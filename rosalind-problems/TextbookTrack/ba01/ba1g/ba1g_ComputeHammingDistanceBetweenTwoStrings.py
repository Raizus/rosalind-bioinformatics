from BioInfoToolkit.Sequences.StringUtils import hamming_distance
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1g/

We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi. For example, CGAAT and CGGAC have two mismatches. The number of mismatches between strings p and q is called the Hamming distance between these strings and is denoted HammingDistance(p, q).

Hamming Distance Problem
    Compute the Hamming distance between two DNA strings.

        Given: Two DNA strings.

        Return: An integer value representing the Hamming distance.
"""

OutputT = int

def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    dist = hamming_distance(seq1, seq2)
    return dist


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    dist = int(lines[0])
    return dist


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    seq1 = lines[0]
    seq2 = lines[1]

    result = solve(seq1, seq2)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1g_1.txt'

    lines = readTextFile(path)
    seq1 = lines[0]
    seq2 = lines[1]

    dist = solve(seq1, seq2)

    out = str(dist)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

