from BioInfoToolkit.Alignment.Alignment import BLOSUM62, SimilarityScore, affineGapGlobalAlignmentCost,globalAlignmentAffineGapPenalty
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/ba5j/

A gap is a contiguous sequence of spaces in a row of an alignment. One way to score gaps more appropriately is to define an affine penalty for a gap of length k as σ + ε · (k - 1), where σ is the gap opening penalty, assessed to the first symbol in the gap, and ε is the gap extension penalty, assessed to each additional symbol in the gap. We typically select ε to be smaller than σ so that the affine penalty for a gap of length k is smaller than the penalty for k independent single-nucleotide indels (σ · k).

Alignment with Affine Gap Penalties Problem
    Construct a highest-scoring global alignment of two strings (with affine gap penalties).

        Given: Two amino acid strings v and w (each of length at most 100).

        Return: The maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score. Use the BLOSUM62 scoring matrix, a gap opening penalty of 11, and a gap extension penalty of 1.
"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res, aligned_seq1_res, aligned_seq2_res = result
    score_sol, _, _ = solution

    is1 = strip_inserts(aligned_seq1_res) == seq1
    is2 = strip_inserts(aligned_seq2_res) == seq2

    openingPenalty: int = -11
    extendingPenalty: int = -1
    similarityScore = SimilarityScore(similarityDict=BLOSUM62)

    score = affineGapGlobalAlignmentCost(
        aligned_seq1_res, aligned_seq2_res, similarityScore, openingPenalty, extendingPenalty)

    correct = is1 == is2 and (score_res == score_sol == score)

    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    openingPenalty: int = -11
    extendingPenalty: int = -1
    similarityScore = SimilarityScore(similarityDict=BLOSUM62)

    # TODO: fix globalAlignmentAffineGapPenalty
    aligned_seq1, aligned_seq2, max_score, _, _, _, _ = globalAlignmentAffineGapPenalty(
        seq1, seq2, openingPenalty, extendingPenalty, similarityScore)

    return max_score, aligned_seq1, aligned_seq2


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    score = int(lines[0])
    aligned_seq1 = lines[1]
    aligned_seq2 = lines[2]

    return score, aligned_seq1, aligned_seq2


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq1 = lines[0]
    seq2 = lines[1]

    result = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, seq1, seq2)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba5j_1.txt'

    lines = readTextFile(path)
    seq1 = lines[0]
    seq2 = lines[1]

    score, aligned_seq1, aligned_seq2 = solve(seq1, seq2)

    print(score)
    print(aligned_seq1)
    print(aligned_seq2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')
    writeTextFile(result_path, aligned_seq1, 'a')
    writeTextFile(result_path, aligned_seq2, 'a')

    correct = solve_and_check(path)
    print(correct)

