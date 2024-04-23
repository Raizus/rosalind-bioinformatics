from BioInfoToolkit.Alignment.Alignment import BLOSUM62, SimilarityScore, globalAlignmentLinearGapPenalty, globalAlignmentLinearGapPenaltyBacktrack, linearGapGlobalAlignmentCost

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/ba5e/

Global Alignment Problem
    Find the highest-scoring alignment between two strings using a scoring matrix.

        Given: Two amino acid strings.

        Return: The maximum alignment score of these strings followed by an alignment achieving this maximum score. Use the BLOSUM62 scoring matrix and indel penalty σ = 5. (If multiple alignments achieving the maximum score exist, you may return any one.)

"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res, aligned_seq1_res, aligned_seq2_res = result
    score_sol, _, _ = solution

    is1 = strip_inserts(aligned_seq1_res) == seq1
    is2 = strip_inserts(aligned_seq2_res) == seq2

    gapPenalty = -5
    similarityScore = SimilarityScore(1, -1, BLOSUM62)

    score = linearGapGlobalAlignmentCost(
        aligned_seq1_res, aligned_seq2_res, gapPenalty, similarityScore)

    correct = is1 == is2 and (score_res == score_sol == score)

    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty = -5
    similarityScore = SimilarityScore(1, -1, BLOSUM62)
    F, _ = globalAlignmentLinearGapPenalty(
        seq1, seq2, gapPenalty, similarityScore)
    aligned_seq1, aligned_seq2 = globalAlignmentLinearGapPenaltyBacktrack(
        seq1, seq2, F, gapPenalty, similarityScore)

    max_score = F[len(seq1)][len(seq2)]

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
    path = f'{cwd}/rosalind_ba5e_1.txt'

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
