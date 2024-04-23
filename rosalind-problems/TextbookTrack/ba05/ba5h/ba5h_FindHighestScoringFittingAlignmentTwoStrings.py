from BioInfoToolkit.Alignment.Alignment import SimilarityScore, linearGapGlobalAlignmentCost, fittingAlignmentLinearGap, fittingAlignmentLinearGapBacktrack

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/ba5h/

Say that we wish to compare the approximately 20,000 amino acid-long NRP synthetase from Bacillus brevis with the approximately 600 amino acid-long A-domain from Streptomyces roseosporus, the bacterium that produces the powerful antibiotic Daptomycin. We hope to find a region within the longer protein sequence v that has high similarity with all of the shorter sequence w. Global alignment will not work because it tries to align all of v to all of w; local alignment will not work because it tries to align substrings of both v and w. Thus, we have a distinct alignment application called the Fitting Alignment Problem.

“Fitting” w to v requires finding a substring v' of v that maximizes the global alignment score between v' and w among all substrings of v.

Fitting Alignment Problem
    Construct a highest-scoring fitting alignment between two strings.

    Given: Two DNA strings v and w, where v has length at most 10000 and w has length at most 1000.

    Return: The maximum score of a fitting alignment of v and w, followed by a fitting alignment achieving this maximum score. Use the simple scoring method in which matches count +1 and both the mismatch and indel penalties are equal to 1. (If multiple fitting alignments achieving the maximum score exist, you may return any one.)
"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res, aligned_seq1_res, aligned_seq2_res = result
    score_sol, _, _ = solution

    is1 = is_subsequence_str(seq1, strip_inserts(aligned_seq1_res))
    is2 = seq2 == strip_inserts(aligned_seq2_res)

    gapPenalty: int = -1
    similarityScore = SimilarityScore(1, -1)
    score = linearGapGlobalAlignmentCost(
        aligned_seq1_res, aligned_seq2_res, gapPenalty, similarityScore)

    correct = is1 == is2 and (score_res == score_sol == score)

    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty: int = -1
    similarityScore = SimilarityScore(1, -1)
    H, _, max_score = fittingAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)
    aligned_seq1, aligned_seq2 = fittingAlignmentLinearGapBacktrack(
        seq1, seq2, H, gapPenalty, similarityScore)

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
    path = f'{cwd}/rosalind_ba5h_1.txt'

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
