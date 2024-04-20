from BioInfoToolkit.Alignment.Alignment import affineGapGlobalAlignmentCost
import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/gaff/

An affine gap penalty is written as a+b⋅(L-1), where L is the length of the gap, a is a positive constant called the gap opening penalty, and b is a positive constant called the gap extension penalty.

We can view the gap opening penalty as charging for the first gap symbol, and the gap extension penalty as charging for each subsequent symbol added to the gap. For example, if a=11 and b=1, then a gap of length 1 would be penalized by 11 (for an average cost of 11 per gap symbol), whereas a gap of length 100 would have a score of 110 (for an average cost of 1.10 per gap symbol).

Consider the strings "PRTEINS" and "PRTWPSEIN". If we use the BLOSUM62 scoring matrix and an affine gap penalty with a=11
and b=1, then we obtain the following optimal alignment.

    PRT---EINS
    |||   |||
    PRTWPSEIN-

Matched symbols contribute a total of 32 to the calculation of the alignment's score, and the gaps cost 13 and 11 respectively, yielding a total score of 8.

    Given: Two protein strings s and t in FASTA format (each of length at most 100 aa).

    Return: The maximum alignment score between s and t, followed by two augmented strings s' and t' representing an optimal alignment of s and t. Use:
        - The BLOSUM62 scoring matrix.
        - Gap opening penalty equal to 11.
        - Gap extension penalty equal to 1.


"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res = result[0]
    score_sol = solution[0]

    res_seq1, res_seq2 = result[1], result[2]
    is1 = is_subsequence_str(seq1, strip_inserts(res_seq1))
    is2 = is_subsequence_str(seq2, strip_inserts(res_seq2))

    openingPenalty: int = -11
    extendingPenalty: int = -1
    similarityScore = Alignment.SimilarityScore(
        similarityDict=Alignment.BLOSUM62)
    cost = affineGapGlobalAlignmentCost(
        alignedSeq1, alignedSeq2, similarityScore, openingPenalty, extendingPenalty)

    correct = (score_res == score_sol == cost) and is1 and is2
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    openingPenalty: int = -11
    extendingPenalty: int = -1

    similarityScore = Alignment.SimilarityScore(
        similarityDict=Alignment.BLOSUM62)
    alignedSeq1, alignedSeq2, score, _, _, _, _ = Alignment.globalAlignmentAffineGapPenalty(
        seq1, seq2, openingPenalty, extendingPenalty, similarityScore)

    return score, alignedSeq1, alignedSeq2


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    score = int(lines[0])
    seq1 = lines[1]
    seq2 = lines[2]

    return score, seq1, seq2


def solve_and_check(input_path: str) -> bool:
    FASTAdict = read_FASTA(input_path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    result = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, seq1, seq2)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_gaff_1.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    score, alignedSeq1, alignedSeq2 = solve(seq1, seq2)

    print(score)
    print(alignedSeq1)
    print(alignedSeq2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')
    writeTextFile(result_path, alignedSeq1, 'a')
    writeTextFile(result_path, alignedSeq2, 'a')

    correct = solve_and_check(path)
    print(correct)

