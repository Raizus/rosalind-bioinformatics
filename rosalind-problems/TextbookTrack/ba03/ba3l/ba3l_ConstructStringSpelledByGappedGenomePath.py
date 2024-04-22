from BioInfoToolkit.Sequences.GenomeAssembly import PairedDeBruijnMultiGraph, deBruijnMultiGraphFromPairedKmers
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba3l/

Gapped Genome Path String Problem
    Reconstruct a string from a sequence of (k,d)-mers corresponding to a path in a paired de Bruijn graph.

        Given: A sequence of (k, d)-mers (a1|b1), ... , (an|bn) such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for all i from 1 to n-1.

        Return: A string Text where the i-th k-mer in Text is equal to Suffix(ai|bi) for all i from 1 to n, if such a string exists.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution

    return correct


def solve(pairs: list[tuple[str, str]], k: int, d: int) -> OutputT:
    _g = deBruijnMultiGraphFromPairedKmers(pairs)
    g = PairedDeBruijnMultiGraph(_g, k, d+k)
    # dot = g.draw_dot()
    # dot.view()

    string = g.reconstructStringFromEulerianPath()
    return string


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    sequence = lines[0]
    return sequence


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, d = [int(v) for v in lines[0].split()]
    pairs = [tuple(a.split('|')[:2]) for a in lines[1:]]
    pairs = [(pair[0], pair[1]) for pair in pairs]

    result = solve(pairs, k, d)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3l_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, d = [int(v) for v in lines[0].split()]
    pairs = [tuple(a.split('|')[:2]) for a in lines[1:]]
    pairs = [(pair[0], pair[1]) for pair in pairs]

    string = solve(pairs, k, d)

    out = string
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
