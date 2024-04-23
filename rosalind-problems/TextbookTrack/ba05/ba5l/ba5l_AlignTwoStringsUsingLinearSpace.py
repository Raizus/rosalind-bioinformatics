from BioInfoToolkit.Alignment.Alignment import BLOSUM62, SimilarityScore, globalAlignmentLinearGapPenaltyScoreInLinearSpace, linearGapGlobalAlignmentCost
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/ba5l/

The pseudocode below for LinearSpaceAlignment describes how to recursively find a longest path in the alignment graph constructed for a substring vtop+1 ... vbottom of v and a substring wleft+1 ... wright of w. LinearSpaceAlignment calls the function MiddleNode(top, bottom, left, right), which returns the coordinate i of the middle node (i, j) defined by the sequences vtop+1 ... vbottom and wleft+1 ... wright. LinearSpaceAlignment also calls MiddleEdge(top, bottom, left, right), which returns → , ↓, or ↘ depending on whether the middle edge is horizontal, vertical, or diagonal. The linear-space alignment of strings v and w is constructed by calling LinearSpaceAlignment(0, n, 0, m). The case left = right describes the alignment of an empty string against the string vtop+1 ... vbottom, which is trivially computed as the score of a gap formed by bottom - top vertical edges.

    LinearSpaceAlignment(top, bottom, left, right)
        if left = right
            return alignment formed by bottom - top vertical edges
        if top = bottom
            return alignment formed by right - left horizontal edges
        middle ← ⌊(left + right)/2⌋

        midNode ← MiddleNode(top, bottom, left, right)
        midEdge ← MiddleEdge(top, bottom, left, right)
        LinearSpaceAlignment(top, midNode, left, middle)
        output midEdge
        if midEdge = "→" or midEdge = "↘"
            middle ← middle + 1
        if midEdge = "↓" or midEdge ="↘"
            midNode ← midNode + 1
        LinearSpaceAlignment(midNode, bottom, middle, right)

Global Alignment in Linear Space Problem
    Find the highest-scoring alignment between two strings using a scoring matrix in linear space.

        Given: Two long amino acid strings (of length approximately 10,000).

        Return: The maximum alignment score of these strings, followed by an alignment achieving this maximum score. Use the BLOSUM62 scoring matrix and indel penalty σ = 5.
"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res, aligned_seq1_res, aligned_seq2_res = result
    score_sol, _, _ = solution

    aux1 = strip_inserts(aligned_seq1_res)
    aux2 = strip_inserts(aligned_seq2_res)
    is1 = aux1 == seq1
    is2 = aux2 == seq2

    gapPenalty = -5
    similarityScore = SimilarityScore(similarityDict=BLOSUM62)
    score = linearGapGlobalAlignmentCost(
        aligned_seq1_res, aligned_seq2_res, gapPenalty, similarityScore)

    correct = is1 == is2 and (score_res == score_sol == score)

    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty = -5
    similarityScore = SimilarityScore(similarityDict=BLOSUM62)
    max_score, aligned_seq1, aligned_seq2 = globalAlignmentLinearGapPenaltyScoreInLinearSpace(
        seq1, seq2, gapPenalty, similarityScore)

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
    path = f'{cwd}/rosalind_ba5l_1.txt'

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
