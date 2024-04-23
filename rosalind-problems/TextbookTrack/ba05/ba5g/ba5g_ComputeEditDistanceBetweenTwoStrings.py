from BioInfoToolkit.Alignment.Alignment import editDistance
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba5g/

In 1966, Vladimir Levenshtein introduced the notion of the edit distance between two strings as the minimum number of edit operations needed to transform one string into another. Here, an edit operation is the insertion, deletion, or substitution of a single symbol. For example, TGCATAT can be transformed into ATCCGAT with five edit operations, implying that the edit distance between these strings is at most 5.

Edit Distance Problem
    Find the edit distance between two strings.

        Given: Two amino acid strings.

        Return: The edit distance between these strings.
"""

OutputT = int


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    dist = editDistance(seq1, seq2)
    return dist


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    dist = int(lines[0])
    return dist


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq1 = lines[0]
    seq2 = lines[1]

    result = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba5g_1.txt'

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

