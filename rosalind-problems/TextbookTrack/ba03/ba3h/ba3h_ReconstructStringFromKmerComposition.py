from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromReads
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba3h/

String Reconstruction Problem
    Reconstruct a string from its k-mer composition.

        Given: An integer k followed by a list of k-mers Patterns.

        Return: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may return any one.)   
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(kmers: list[str], k: int) -> OutputT:
    _g = deBruijnMultiGraphFromReads(kmers, k)
    g = DeBruijnMultiGraph(_g, k)
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
    k = int(lines[0])
    kmers = lines[1:]

    result = solve(kmers, k)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3h_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k = int(lines[0])
    kmers = lines[1:]

    string = solve(kmers, k)

    out = string
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
