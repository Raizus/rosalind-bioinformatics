from BioInfoToolkit.Alignment.Alignment import globalAlignmentLinearGapPenaltyBacktrack, linearGapGlobalAlignmentCost

import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/glob/

To penalize symbol substitutions differently depending on which two symbols are involved in the substitution, we obtain a scoring matrix S in which Si,j represents the (negative) score assigned to a substitution of the ith symbol of our alphabet A with the jth symbol of A.

A gap penalty is the component deducted from alignment score due to the presence of a gap. A gap penalty may be a function of the length of the gap; for example, a linear gap penalty is a constant g such that each inserted or deleted symbol is charged g; as a result, the cost of a gap of length L is equal to gL.

    Given: Two protein strings s and t in FASTA format (each of length at most 1000 aa).

    Return: The maximum alignment score between s and t. Use:

    - The BLOSUM62 scoring matrix.
    - Linear gap penalty equal to 5 (i.e., a cost of -5 is assessed for each gap symbol).

"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> int:
    gapPenalty = -5
    similarityScore = Alignment.SimilarityScore(1, -1, Alignment.BLOSUM62)
    F, _ = Alignment.globalAlignmentLinearGapPenalty(
        seq1, seq2, gapPenalty, similarityScore)
    
    max_score = F[len(seq1)][len(seq2)]
    return max_score


def load_results(path: str) -> int:
    lines = readTextFile(path)
    score = int(lines[0])

    return score


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
    path = f'{cwd}/rosalind_glob_1.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    gapPenalty = -5
    similarityScore = Alignment.SimilarityScore(1, -1, Alignment.BLOSUM62)
    F, _ = Alignment.globalAlignmentLinearGapPenalty(
        seq1, seq2, gapPenalty, similarityScore)
    alignedSeq1, alignedSeq2 = globalAlignmentLinearGapPenaltyBacktrack(
        seq1, seq2, F, gapPenalty, similarityScore)

    max_score = F[len(seq1)][len(seq2)]
    print(max_score)

    # score2 = linearGapGlobalAlignmentCost(
    #     alignedSeq1, alignedSeq2, gapPenalty, similarityScore)
    # print(score2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(max_score), 'w')

    correct = solve_and_check(path)
    print(correct)
