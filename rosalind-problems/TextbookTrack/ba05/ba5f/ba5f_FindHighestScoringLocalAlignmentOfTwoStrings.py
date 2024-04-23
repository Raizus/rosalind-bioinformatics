from BioInfoToolkit.Alignment.Alignment import PAM250, SimilarityScore, linearGapGlobalAlignmentCost, localAlignmentLinearGap

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/ba5f/

Local Alignment Problem
    Find the highest-scoring local alignment between two strings.

    Given: Two amino acid strings.

    Return: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum score. Use the PAM250 scoring matrix and indel penalty σ = 5. (If multiple local alignments achieving the maximum score exist, you may return any one.)

"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res, aligned_seq1_res, aligned_seq2_res = result
    score_sol, _, _ = solution

    is1 = is_subsequence_str(seq1, strip_inserts(aligned_seq1_res))
    is2 = is_subsequence_str(seq2, strip_inserts(aligned_seq2_res))

    gapPenalty = -5
    similarityScore = SimilarityScore(1, -1, PAM250)
    score = linearGapGlobalAlignmentCost(
        aligned_seq1_res, aligned_seq2_res, gapPenalty, similarityScore)

    correct = is1 == is2 and (score_res == score_sol == score)

    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty = -5
    similarityScore = SimilarityScore(1, -1, PAM250)
    aligned_seq1, aligned_seq2, H, max_ind = localAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)
    i, j = max_ind[0]
    max_score = H[i][j]
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
    path = f'{cwd}/rosalind_ba5f_1.txt'

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

