from BioInfoToolkit.Sequences.BarrowsWheeler import lastToFirst, lastToFirstArray
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9k/

The Last-to-First array, denoted LastToFirst(i), answers the following question: given a symbol at position i in LastColumn, what is its position in FirstColumn?

Last-to-First Mapping Problem

    Given: A string Transform and an integer i.

    Return: The position LastToFirst(i) in FirstColumn in the Burrows-Wheeler matrix if LastColumn = Transform.
"""

OutputT = int


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(bwt: str, i: int) -> OutputT:
    # last_to_first_array = lastToFirstArray(bwt)
    # k = last_to_first_array[i]
    k = lastToFirst(bwt, i)
    return k


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    pos = int(lines[0])
    return pos


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    bwt = lines[0]
    i = int(lines[1])

    result = solve(bwt, i)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9k_1.txt'

    lines = readTextFile(path)
    bwt = lines[0]
    i = int(lines[1])

    pos = solve(bwt, i)

    out = str(pos)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
