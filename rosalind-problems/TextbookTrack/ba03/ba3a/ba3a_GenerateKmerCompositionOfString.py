from BioInfoToolkit.Sequences.StringUtils import kmer_composition
from collections import Counter
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba3a/

Given a string Text, its k-mer composition Compositionk(Text) is the collection of all k-mer substrings of Text (including repeated k-mers). For example,

    Composition3(TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}

Note that we have listed k-mers in lexicographic order (i.e., how they would appear in a dictionary) rather than in the order of their appearance in TATGGGGTGC. We have done this because the correct ordering of the reads is unknown when they are generated.

String Composition Problem
    Generate the k-mer composition of a string.

        Given: An integer k and a string Text.

        Return: Compositionk(Text) (the k-mers can be provided in any order).
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(sequence: str, k: int) -> OutputT:
    kmers = kmer_composition(sequence, k)
    return kmers


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]
    return kmers


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k = int(lines[0])
    sequence = lines[1]

    result = solve(sequence, k)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3a_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k = int(lines[0])
    sequence = lines[1]

    kmers = solve(sequence, k)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for kmer in kmers:
        print(kmer)
        writeTextFile(result_path, kmer, 'a')

    correct = solve_and_check(path)
    print(correct)
