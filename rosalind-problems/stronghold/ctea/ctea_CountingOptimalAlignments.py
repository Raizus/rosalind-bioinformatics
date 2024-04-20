import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.Alignment.Alignment import editDistance, hamming_distance
from BioInfoToolkit.IO import read_FASTA

import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ctea/

Recall from “Edit Distance Alignment” that if s' and t' are the augmented strings corresponding to an alignment of strings s and t, then the edit alignment score of s' and t' was given by the Hamming distance dH(s',t') (because s' and t'have the same length and already include gap symbols to denote insertions/deletions).

As a result, we obtain dE(s,t)=mins',t'dH(s',t'), where the minimum is taken over all alignments of s and t. Strings s' and t' achieving this minimum correspond to an optimal alignment with respect to edit alignment score.

    Given: Two protein strings sand t in FASTA format, each of length at most 1000 aa.

    Return: The total number of optimal alignments of s and t with respect to edit alignment score, modulo 134,217,727 (227-1).
    """


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> int:
    mod = 2**27 - 1
    gapPenalty = -1
    similarityScore = Alignment.SimilarityScore(0, -1)
    count = Alignment.countOptimalGlobalAlignmentsLinearGap(
        seq1, seq2, gapPenalty, similarityScore, mod)
    return count


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
    path = f'{cwd}/rosalind_ctea_1.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    mod = 2**27 - 1
    gapPenalty = -1
    similarityScore = Alignment.SimilarityScore(0, -1)
    count = Alignment.countOptimalGlobalAlignmentsLinearGap(
        seq1, seq2, gapPenalty, similarityScore, mod)

    # dist = Alignment.editDistance(seq1, seq2)
    # print(dist)
    # F, _ = Alignment.globalAlignmentLinearGapPenalty(
    #     seq1, seq2, gapPenalty, similarityScore)
    # alignedSeq1, alignedSeq2 = Alignment.globalAlignmentLinearGapPenaltyBacktrack(
    #     seq1, seq2, F, gapPenalty, similarityScore)
    # hamm_dist = hamming_distance(alignedSeq1, alignedSeq2)
    # print("Hamming Distance: ", hamm_dist)
    # print("Edit distance", editDistance(alignedSeq1, alignedSeq2))


    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    correct = solve_and_check(path)
    print(correct)
