from BioInfoToolkit.Sequences.SequenceUtils import reverseComplement
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1c/

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'. Given a nucleotide p, we denote its complementary nucleotide as p. The reverse complement of a DNA string Pattern = p1…pn is the string Pattern = pn … p1 formed by taking the complement of each nucleotide in Pattern, then reversing the resulting string.

For example, the reverse complement of Pattern = "GTCA" is Pattern = "TGAC". 

Reverse Complement Problem
    Find the reverse complement of a DNA string.

        Given: A DNA string Pattern.

        Return: Pattern, the reverse complement of Pattern.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(string: str) -> OutputT:
    rcomp = reverseComplement(string, 'DNA')
    return rcomp


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    rcomp = lines[0]
    return rcomp


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    string = lines[0]
    result = solve(string)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1c_1.txt'

    lines = readTextFile(path)
    string = lines[0]

    rcomp = solve(string)

    print(rcomp)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, rcomp, 'w')

    correct = solve_and_check(path)
    print(correct)
