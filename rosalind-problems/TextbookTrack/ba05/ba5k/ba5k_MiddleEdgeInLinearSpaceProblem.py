from BioInfoToolkit.Alignment.Alignment import BLOSUM62, SimilarityScore, middle_edge
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import re

"""
https://rosalind.info/problems/ba5k/

Middle Edge in Linear Space Problem
    Find a middle edge in the alignment graph in linear space.

        Given: Two amino acid strings.

        Return: A middle edge in the alignment graph of these strings, where the optimal path is defined by the BLOSUM62 scoring matrix and a linear indel penalty equal to 5. Return the middle edge in the form “(i, j) (k, l)”, where (i, j) connects to (k, l).
"""

OutputT = tuple[tuple[int, int], tuple[int, int]]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty = -5
    similarityScore = SimilarityScore(similarityDict=BLOSUM62)

    edge, _ = middle_edge(seq1, seq2, gapPenalty, similarityScore)

    return edge


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    match = re.match(r'\((\d+),\s?(\d+)\)\s+\((\d+),\s?(\d+)\)', lines[0])
    if match is None:
        raise Exception(f'{lines[0]} does not have the correct format.')

    edge: OutputT = ((int(match[1]), int(match[2])),
                     (int(match[3]), int(match[4])))

    return edge


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq1 = lines[0]
    seq2 = lines[1]

    result = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba5k_1.txt'

    lines = readTextFile(path)
    seq1 = lines[0]
    seq2 = lines[1]

    edge = solve(seq1, seq2)

    out = f"({edge[0][0]}, {edge[0][1]}) ({edge[1][0]}, {edge[1][1]})"

    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

