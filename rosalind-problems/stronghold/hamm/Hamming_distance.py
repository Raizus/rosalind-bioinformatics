
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import hamming_distance

"""
Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), is the number of corresponding symbols that differ in s and t.

    Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

    Return: The Hamming distance d_H(s,t).
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> int:
    dist = hamming_distance(seq1, seq2)
    return dist


def load_results(path: str) -> int:
    lines = readTextFile(path)
    dist = int(lines[0])
    return dist


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    seq1 = lines[0]
    seq2 = lines[1]

    dist = solve(seq1, seq2)
    solution = load_results(solution_path)

    correct = verify(dist, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_hamm_1.txt'

    lines = readTextFile(path)
    seq1 = lines[0]
    seq2 = lines[1]

    dist = solve(seq1, seq2)
    print(dist)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(dist), 'w')


