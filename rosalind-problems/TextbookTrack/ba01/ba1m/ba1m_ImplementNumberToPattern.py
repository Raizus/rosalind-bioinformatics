from BioInfoToolkit.Sequences.StringUtils import numberToPattern
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1m/

Convert an integer to its corresponding DNA string.

    Given: Integers index and k.

    Return: NumberToPattern(index, k).
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(index: int, k: int) -> OutputT:
    alphabet = ['A', 'C', 'G', 'T']
    pattern = numberToPattern(index, k, alphabet)
    return pattern


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    pattern = lines[0]
    return pattern


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    index = int(lines[0])
    k = int(lines[1])

    result = solve(index, k)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1m_1.txt'

    lines = readTextFile(path)
    index = int(lines[0])
    k = int(lines[1])

    pattern = solve(index, k)

    out = pattern
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
