from collections import Counter
from BioInfoToolkit.Sequences.StringUtils import generate_d_neighborhood
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1n/

The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose Hamming distance from Pattern does not exceed d.
Generate the d-Neighborhood of a String

Find all the neighbors of a pattern.

    Given: A DNA string Pattern and an integer d.

    Return: The collection of strings Neighbors(Pattern, d).
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(pattern: str, d: int) -> OutputT:
    alphabet = ['A', 'C', 'G', 'T']
    kmers: OutputT = []
    for kmer in generate_d_neighborhood(pattern, d, set(alphabet)):
        kmers.append(kmer)
    return kmers


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]
    return kmers


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    pattern = lines[0]
    d = int(lines[1])

    result = solve(pattern, d)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1n_1.txt'

    lines = readTextFile(path)
    pattern = lines[0]
    d = int(lines[1])

    strings = solve(pattern, d)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for kmer in strings:
        print(kmer)
        writeTextFile(result_path, kmer, 'a')

    correct = solve_and_check(path)
    print(correct)
