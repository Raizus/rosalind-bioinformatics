from BioInfoToolkit.Sequences.SequenceUtils import minimum_skew
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1f/

Define the skew of a DNA string Genome, denoted Skew(Genome), as the difference between the total number of occurrences of 'G' and 'C' in Genome. Let Prefixi (Genome) denote the prefix (i.e., initial substring) of Genome of length i. For example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCC")) are:

    0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem
    Find a position in a genome minimizing the skew.

        Given: A DNA string Genome.

        Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i (from 0 to |Genome|).
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(sequence: str) -> OutputT:
    _, positions = minimum_skew(sequence)
    return positions


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    idxs = [int(v) for v in lines[0].split()]
    return idxs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    sequence = lines[0]
    result = solve(sequence)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1f_1.txt'

    lines = readTextFile(path)
    sequence = lines[0]

    positions = solve(sequence)

    out = ' '.join(str(p) for p in positions)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
