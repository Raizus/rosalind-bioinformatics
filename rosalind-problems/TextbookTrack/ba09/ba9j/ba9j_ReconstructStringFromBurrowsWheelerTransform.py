from BioInfoToolkit.Sequences.BarrowsWheeler import stringFromBWT
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9j/

In “Construct the Burrows-Wheeler Transform of a String”, we introduced the Burrows-Wheeler transform of a string Text. In this problem, we give you the opportunity to reverse this transform.

Inverse Burrows-Wheeler Transform Problem
    Reconstruct a string from its Burrows-Wheeler transform.

        Given: A string Transform (with a single "$" sign).

        Return: The string Text such that BWT(Text) = Transform.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(bwt: str) -> OutputT:
    string = stringFromBWT(bwt)
    return string


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    string = lines[0]
    return string


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    bwt = lines[0]

    result = solve(bwt)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9j_1.txt'

    lines = readTextFile(path)
    bwt = lines[0]

    string = solve(bwt)

    out = string
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
