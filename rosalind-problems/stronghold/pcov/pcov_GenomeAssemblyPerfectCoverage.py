from collections import Counter
from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromReads
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import kmer_gen

"""
    https://rosalind.info/problems/pcov/

A circular string is a string that does not have an initial or terminal element; instead, the string is viewed as a necklace of symbols. We can represent a circular string as a string enclosed in parentheses. For example, consider the circular DNA string (ACGTAC), and note that because the string "wraps around" at the end, this circular string can equally be represented by (CGTACA), (GTACAC), (TACACG), (ACACGT), and (CACGTA). The definitions of substrings and superstrings are easy to generalize to the case of circular strings (keeping in mind that substrings are allowed to wrap around).

    Given: A collection of (error-free) DNA k-mers (k≤50) taken from the same strand of a circular chromosome. In this dataset, all k-mers from this strand of the chromosome are present, and their de Bruijn graph consists of exactly one simple cycle.

    Return: A cyclic superstring of minimal length containing the reads (thus corresponding to a candidate cyclic chromosome).
"""


def verify(result: str, solution: str, kmers: list[str]) -> bool:
    k = len(kmers[0])
    result_kmers = Counter(kmer for kmer in kmer_gen(result, k, True))

    correct = len(result) == len(solution) and result_kmers == Counter(kmers)
    return correct


def solve(kmers: list[str]) -> str:
    k = len(kmers[0])

    _g = deBruijnMultiGraphFromReads(kmers)
    g = DeBruijnMultiGraph(_g, k)

    string = g.reconstructStringFromEulerianCycle()
    return string


def load_results(path: str) -> str:
    lines = readTextFile(path)
    string = lines[0]
    return string


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    kmers = [line for line in lines if len(line) and not line.isspace()]
    string = solve(kmers)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(string, solution, kmers)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_pcov_1.txt'

    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]
    k = len(kmers[0])

    _g = deBruijnMultiGraphFromReads(kmers)
    g = DeBruijnMultiGraph(_g, k)
    # dot = g.draw_dot()
    # dot.view()

    string = g.reconstructStringFromEulerianCycle()
    print(string)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, string, 'w')

    correct = solve_and_check(path)
    print(correct)
