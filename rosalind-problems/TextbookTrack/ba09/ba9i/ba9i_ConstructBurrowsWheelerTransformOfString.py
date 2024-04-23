from BioInfoToolkit.Sequences.BarrowsWheeler import barrowsWheelerTransform
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9i/

Our goal is to further improve on the memory requirements of the suffix array introduced in “Construct the Suffix Array of a String” for multiple pattern matching.

Given a string Text, form all possible cyclic rotations of Text; a cyclic rotation is defined by chopping off a suffix from the end of Text and appending this suffix to the beginning of Text. Next — similarly to suffix arrays — order all the cyclic rotations of Text lexicographically to form a |Text| x |Text| matrix of symbols that we call the Burrows-Wheeler matrix and denote by M(Text).

Note that the first column of M(Text) contains the symbols of Text ordered lexicographically. In turn, the second column of M(Text) contains the second symbols of all cyclic rotations of Text, and so it too represents a (different) rearrangement of symbols from Text. The same reasoning applies to show that any column of M(Text) is some rearrangement of the symbols of Text. We are interested in the last column of M(Text), called the Burrows-Wheeler transform of Text, or BWT(Text).

Burrows-Wheeler Transform Construction Problem
    Construct the Burrows-Wheeler transform of a string.

        Given: A string Text.

        Return: BWT(Text).
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(text: str) -> OutputT:
    bwt = barrowsWheelerTransform(text)
    return bwt


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    bwt = lines[0]
    return bwt


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    string = lines[0]

    result = solve(string)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9i_1.txt'

    lines = readTextFile(path)
    text = lines[0]

    bwt = solve(text)

    out = bwt
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

