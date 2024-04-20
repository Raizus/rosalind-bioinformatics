from BioInfoToolkit.Sequences.StringUtils import getMedianStrings, hamming_distance
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba2b/

Given a k-mer Pattern and a longer string Text, we use d(Pattern, Text) to denote the minimum Hamming distance between Pattern and any k-mer in Text,

    d(Pattern,Text)=min all k-mers Pattern' in TextHammingDistance(Pattern,Pattern').

Given a k-mer Pattern and a set of strings Dna = {Dna1, … , Dnat}, we define d(Pattern, Dna) as the sum of distances between Pattern and all strings in Dna,

    d(Pattern,Dna)=∑i=1td(Pattern,Dnai).

Our goal is to find a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern, the same task that the Equivalent Motif Finding Problem is trying to achieve. We call such a k-mer a median string for Dna.

Median String Problem
    Find a median string.

        Given: An integer k and a collection of strings Dna.

        Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern. (If multiple answers exist, you may return any one.)
"""

OutputT = str


def verify(result: OutputT, solution: OutputT, sequences: list[str]) -> bool:
    dist_result = sum([hamming_distance(result, seq) for seq in sequences])
    dist_solution = sum([hamming_distance(solution, seq) for seq in sequences])
    correct = len(result) == len(solution) and dist_result == dist_solution
    return correct


def solve(sequences: list[str], k: int) -> OutputT:
    alphabet = set(['A', 'C', 'G', 'T'])
    median_strings, _ = getMedianStrings(sequences, k, alphabet)
    median_string = next(iter(median_strings))
    return median_string


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    median_string = lines[0]
    return median_string


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k = int(lines[0])
    sequences = lines[1:]

    result = solve(sequences, k)

    solution = load_results(solution_path)

    correct = verify(result, solution, sequences)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba2b_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k = int(lines[0])
    sequences = lines[1:]

    string = solve(sequences, k)

    out = string
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
