from collections import Counter
from BioInfoToolkit.Sequences.StringUtils import generate_kmers_from_alphabet, kmer_gen
from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromReads
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba3i/

A k-universal circular string is a circular string that contains every possible k-mer constructed over a given alphabet.

k-Universal Circular String Problem
    Find a k-universal circular binary string.

        Given: An integer k.

        Return: A k-universal circular string. (If multiple answers exist, you may return any one.)
"""

OutputT = str


def verify(result: OutputT, solution: OutputT, k: int) -> bool:
    alphabet = ['0', '1']
    alphabet_kmers = generate_kmers_from_alphabet(k, alphabet)
    string_kmers = kmer_gen(result, k, True)

    correct = Counter(alphabet_kmers) == Counter(string_kmers)

    return correct


def solve(k: int) -> OutputT:
    alphabet = ['0', '1']
    kmers_gen = generate_kmers_from_alphabet(k, alphabet)

    _g = deBruijnMultiGraphFromReads(kmers_gen, k)
    g = DeBruijnMultiGraph(_g, k)
    # dot = g.draw_dot()
    # dot.view()

    string = g.reconstructStringFromEulerianCycle()
    return string


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    sequence = lines[0]
    return sequence


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    k = int(lines[0])

    result = solve(k)

    solution = load_results(solution_path)

    correct = verify(result, solution, k)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3i_1.txt'

    lines = readTextFile(path)
    k = int(lines[0])

    string = solve(k)

    out = string
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
