from BioInfoToolkit.Sequences.StringUtils import partialSuffixArray
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9q/

To construct the partial suffix array SuffixArrayk(Text), we first need to construct the full suffix array and then retain only the elements of this array that are divisible by K, along with their indices i.

Partial Suffix Array Construction Problem
    Construct the partial suffix array of a string.

        Given: A string Text and a positive integer K.

        Return: SuffixArrayK(Text), in the form of a list of ordered pairs (i, SuffixArray(i)) for all nonempty entries in the partial suffix array.
"""

OutputT = list[tuple[int, int]]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(string: str, k: int) -> OutputT:
    partial_array = partialSuffixArray(string, k)
    return partial_array


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    partial_array: OutputT = []
    for line in lines:
        v1, v2 = line.split(',')
        partial_array.append((int(v1), int(v2)))
    return partial_array


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    string = lines[0]
    k = int(lines[1])

    result = solve(string, k)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9q_1.txt'

    lines = readTextFile(path)
    string = lines[0]
    k = int(lines[1])

    partial_array = solve(string, k)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for i, v in partial_array:
        out = f"{i},{v}"
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
