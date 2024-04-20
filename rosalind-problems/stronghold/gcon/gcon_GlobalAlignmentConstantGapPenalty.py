import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/gcon/

In a constant gap penalty, every gap receives some predetermined constant penalty, regardless of its length. Thus, the insertion or deletion of 1000 contiguous symbols is penalized equally to that of a single symbol.

    Given: Two protein strings s and t in FASTA format (each of length at most 1000 aa).

    Return: The maximum alignment score between s and t. Use:
        - The BLOSUM62 scoring matrix.
        - Constant gap penalty equal to 5.

    """


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> int:
    gapPenalty = -5
    similarityScore = Alignment.SimilarityScore(
        similarityDict=Alignment.BLOSUM62)
    score, _ = Alignment.globalAlignmentConstantGapPenalty(
        seq1, seq2, gapPenalty, similarityScore)
    return score


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])

    return count


def solve_and_check(input_path: str) -> bool:
    FASTAdict = read_FASTA(input_path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    result = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_gcon_2.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    score = solve(seq1, seq2)
    print(score)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')

    # correct = solve_and_check(path)
    # print(correct)
