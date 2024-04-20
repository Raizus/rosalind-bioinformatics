from BioInfoToolkit.Sequences.StringUtils import findMotif
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1d/

In this problem, we ask a simple question: how many times can one string occur as a substring of another? Recall from “Find the Most Frequent Words in a String” that different occurrences of a substring can overlap with each other. For example, ATA occurs three times in CGATATATCCATAG.

Pattern Matching Problem
    Find all occurrences of a pattern in a string.

        Given: Strings Pattern and Genome.

        Return: All starting positions in Genome where Pattern appears as a substring. Use 0-based indexing.
"""

OutputT = list[int]

def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(sequence: str, pattern: str) -> OutputT:
    matches = findMotif(pattern, sequence)
    matches = sorted(matches.keys())
    return matches


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    idxs = [int(v) for v in lines[0].split()]
    return idxs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    pattern = lines[0]
    sequence = lines[1]
    result = solve(sequence, pattern)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1d_1.txt'

    lines = readTextFile(path)
    pattern = lines[0]
    sequence = lines[1]

    matches = solve(sequence, pattern)

    out = ' '.join(str(m) for m in matches)

    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
