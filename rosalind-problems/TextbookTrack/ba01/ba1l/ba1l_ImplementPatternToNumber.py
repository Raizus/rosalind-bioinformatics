from BioInfoToolkit.Sequences.StringUtils import patternToNumber
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1l/

Convert a DNA string to a number.

    Given: A DNA string Pattern.

    Return: PatternToNumber(Pattern).
"""

OutputT = int


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(sequence: str) -> OutputT:
    alphabet = ['A', 'C', 'G', 'T']
    numb = patternToNumber(sequence, alphabet)
    return numb


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    numb = int(lines[0])
    return numb


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    sequence = lines[0]

    result = solve(sequence)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1l_0.txt'

    lines = readTextFile(path)
    sequence = lines[0]

    numb = solve(sequence)

    out = str(numb)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
