from BioInfoToolkit.Alignment.Alignment import SimilarityScore, overlap_trim_s_t, overlapAlignmentLinearGap, linearGapGlobalAlignmentCost
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/ba5i/

When we assembled genomes, we discussed how to use overlapping reads to assemble a genome, a problem that was complicated by errors in reads. We would like to find overlaps between error-prone reads as well.

An overlap alignment of strings v = v1 ... vn and w = w1 ... wm is a global alignment of a suffix of v with a prefix of w. An optimal overlap alignment of strings v and w maximizes the global alignment score between an i-suffix of v and a j-prefix of w (i.e., between vi ... vn and w1 ... wj) among all i and j.

Overlap Alignment Problem
    Construct a highest-scoring overlap alignment between two strings.

        Given: Two protein strings v and w, each of length at most 1000.

        Return: The score of an optimal overlap alignment of v and w, followed by an alignment of a suffix v' of v and a prefix w' of w achieving this maximum score. Use an alignment score in which matches count +1 and both the mismatch and indel penalties are 2. (If multiple overlap alignments achieving the maximum score exist, you may return any one.)
"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res, aligned_seq1_res, aligned_seq2_res = result
    score_sol, _, _ = solution

    aux1 = strip_inserts(aligned_seq1_res)
    aux2 = strip_inserts(aligned_seq2_res)
    is1 = aux1 == seq1[-len(aux1):]
    is2 = aux2 == seq2[:len(aux2)]

    gapPenalty: int = -2
    similarityScore = SimilarityScore(1, -2)
    score = linearGapGlobalAlignmentCost(
        aligned_seq1_res, aligned_seq2_res, gapPenalty, similarityScore)

    correct = is1 == is2 and (score_res == score_sol == score)

    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty: int = -2
    similarityScore = SimilarityScore(1, -2)
    aligned_seq1, aligned_seq2, _, max_score = overlapAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)
    
    new_s, new_t = overlap_trim_s_t(aligned_seq1, aligned_seq2)

    return max_score, new_s, new_t


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
    path = f'{cwd}/rosalind_ba5i_1.txt'

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
