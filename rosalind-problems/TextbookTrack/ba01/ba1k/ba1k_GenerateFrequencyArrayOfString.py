from BioInfoToolkit.Sequences.StringUtils import kmers_frequency_array
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1k/

Given an integer k, we define the frequency array of a string Text as an array of length 4k, where the i-th element of the array holds the number of times that the i-th k-mer (in the lexicographic order) appears in Text (see Figure 1.

Computing a Frequency Array
    Generate the frequency array of a DNA string.

        Given: A DNA string Text and an integer k.

        Return: The frequency array of k-mers in Text.
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(sequence: str, k: int) -> OutputT:
    alphabet = ['A', 'C', 'G', 'T']
    array = kmers_frequency_array(sequence, k, alphabet)
    return array


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    array = [int(v) for v in lines[0].split()]
    return array


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    sequence = lines[0]
    k = int(lines[1])

    result = solve(sequence, k)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1k_1.txt'

    lines = readTextFile(path)
    sequence = lines[0]
    k = int(lines[1])

    array = solve(sequence, k)

    out = ' '.join(str(v) for v in array)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
