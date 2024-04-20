from BioInfoToolkit.Sequences.StringUtils import findClumps
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1e/

Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger) string Genome if there is an interval of Genome of length L in which Pattern appears at least t times. For example, TGCA forms a (25,3)-clump in the following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

Clump Finding Problem
    Find patterns forming clumps in a string.

        Given: A string Genome, and integers k, L, and t.

        Return: All distinct k-mers forming (L, t)-clumps in Genome.
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(sequence: str, k: int, L: int, t: int) -> OutputT:
    clumps = findClumps(sequence, k, L, t)
    return list(clumps)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    kmers = lines[0].split()
    return kmers


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    sequence = lines[0]
    k, L, t = [int(v) for v in lines[1].split()]
    result = solve(sequence, k, L, t)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1e_1.txt'

    lines = readTextFile(path)
    sequence = lines[0]
    k, L, t = [int(v) for v in lines[1].split()]
    clumps = solve(sequence, k, L, t)

    out = ' '.join(clumps)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
