import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/edta/

An alignment of two strings s and t is defined by two strings s' and t' satisfying the following three conditions: 1. s' and t' must be formed from adding gap symbols "-" to each of s and t, respectively; as a result, s and t will form subsequences of s' and t'. 2. s' and t' must have the same length. 3. Two gap symbols may not be aligned; that is, if s'[j] is a gap symbol, then t'[j] cannot be a gap symbol, and vice-versa.

We say that s' and t' augment s and t. Writing s' directly over t' so that symbols are aligned provides us with a scenario for transforming s into t. Mismatched symbols from s and t correspond to symbol substitutions; a gap symbol s'[j] aligned with a non-gap symbol t'[j] implies the insertion of this symbol into t; a gap symbol t'[j] aligned with a non-gap symbol s'[j] implies the deletion of this symbol from s.

Thus, an alignment represents a transformation of s into t via edit operations. We define the corresponding edit alignment score of s' and t' as dH(s',t') (Hamming distance is used because the gap symbol has been introduced for insertions and deletions). It follows that dE(s,t)=mins',t'dH(s',t'), where the minimum is taken over all alignments of s and t. We call such a minimum score alignment an optimal alignment (with respect to edit distance).

    Given: Two protein strings s and t in FASTA format (with each string having length at most 1000 aa).

    Return: The edit distance dE(s,t) followed by two augmented strings s' and t' representing an optimal alignment of s and t.
    """


def verify(result: tuple[int, str, str], solution: tuple[int, str, str]) -> bool:
    correct = result[0] == solution[0] == Alignment.editDistance(
        result[1], result[2])
    return correct


def solve(seq1: str, seq2: str) -> tuple[int, str, str]:
    dist = Alignment.editDistance(seq1, seq2)

    matchScore = 0
    similarityScore = Alignment.SimilarityScore(
        matchScore, similarityDict=None)
    gapPenalty = -1
    F, _ = Alignment.globalAlignmentLinearGapPenalty(
        seq1, seq2, gapPenalty, similarityScore)
    alignedSeq1, alignedSeq2 = Alignment.globalAlignmentLinearGapPenaltyBacktrack(
        seq1, seq2, F, gapPenalty, similarityScore)
    return dist, alignedSeq1, alignedSeq2


def load_results(path: str) -> tuple[int, str, str]:
    lines = readTextFile(path)
    dist = int(lines[0])
    seq1, seq2 = lines[1], lines[2]

    return dist, seq1, seq2


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
    path = f'{cwd}/rosalind_edta_1.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    dist, alignedSeq1, alignedSeq2 = solve(seq1, seq2)

    print(dist)
    print(alignedSeq1)
    print(alignedSeq2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(dist), 'w')
    writeTextFile(result_path, alignedSeq1, 'a')
    writeTextFile(result_path, alignedSeq2, 'a')

    correct = solve_and_check(path)
    print(correct)