from BioInfoToolkit.Sequences.StringUtils import suffixArray
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9g/

In “Find the Longest Repeat in a String”, we saw that suffix trees can be too memory intensive to apply in practice.

In 1993, Udi Manber and Gene Myers introduced suffix arrays as a memory-efficient alternative to suffix trees. To construct SuffixArray(Text), we first sort all suffixes of Text lexicographically, assuming that "$" comes first in the alphabet. The suffix array is the list of starting positions of these sorted suffixes. For example,

    SuffixArray("panamabananas$") = (13, 5, 3, 1, 7, 9, 11, 6, 4, 2, 8, 10, 0, 12).

Suffix Array Construction Problem
    Construct the suffix array of a string.

        Given: A string Text.

        Return: SuffixArray(Text).
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(string: str) -> OutputT:
    array = suffixArray(string)
    return array


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    array = [int(v) for v in lines[0].split(',')]
    return array


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
    path = f'{cwd}/rosalind_ba9g_1.txt'

    lines = readTextFile(path)
    string = lines[0]

    array = solve(string)

    out = ', '.join(str(v) for v in array)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
