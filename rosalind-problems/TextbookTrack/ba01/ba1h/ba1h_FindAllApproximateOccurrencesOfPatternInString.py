from BioInfoToolkit.Sequences.StringUtils import approximatePatternMatching
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1h/

We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA box may appear with slight variations leads to the following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem
    Find all approximate occurrences of a pattern in a string.

        Given: Strings Pattern and Text along with an integer d.

        Return: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution) and len(result) == len(solution)
    return correct


def solve(sequence: str, pattern: str, d: int) -> OutputT:
    idxs = approximatePatternMatching(sequence, pattern, d)
    return idxs


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    idxs = [int(v) for v in lines[0].split()]
    return idxs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    pattern = lines[0]
    sequence = lines[1]
    d = int(lines[2])

    result = solve(sequence, pattern, d)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1h_1.txt'

    lines = readTextFile(path)
    pattern = lines[0]
    sequence = lines[1]
    d = int(lines[2])

    idxs = solve(sequence, pattern, d)

    out = ' '.join(str(i) for i in idxs)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
