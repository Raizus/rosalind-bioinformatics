from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba3b/

String Spelled by a Genome Path Problem
    Find the string spelled by a genome path.

        Given: A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1 symbols of Patterni are equal to the first k - 1 symbols of Patterni+1 for i from 1 to n-1.

        Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.
"""

OutputT = str


def stitchKmers(kmers: list[str]) -> str:
    sequence = kmers[0]
    for kmer in kmers[1:]:
        sequence += kmer[-1]
    return sequence


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(kmers: list[str]) -> OutputT:
    sequence = stitchKmers(kmers)
    return sequence


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    sequence = lines[0]
    return sequence


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    result = solve(kmers)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3b_1.txt'

    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    sequence = solve(kmers)

    out = sequence
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
