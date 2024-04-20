import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/loca/

A local alignment of two strings s and t is an alignment of substrings r and u of s and t, respectively. Let opt(r,u) denote the score of an optimal alignment of r and u with respect to some predetermined alignment score.

    Given: Two protein strings s and t in FASTA format (each having length at most 1000 aa).

    Return: A maximum alignment score along with substrings r and u of s and t, respectively, which produce this maximum alignment score (multiple solutions may exist, in which case you may output any one). Use:
        - The PAM250 scoring matrix.
        - Linear gap penalty equal to 5.

"""


def verify(result: tuple[int, str, str], solution: tuple[int, str, str], seq1: str, seq2: str) -> bool:
    score_res = result[0]
    score_sol = solution[0]

    res_seq1, res_seq2 = result[1], result[2]
    is1 = is_subsequence_str(seq1, res_seq1)
    is2 = is_subsequence_str(seq2, res_seq2)

    # should confirm score aligned sequences, but insertions were stripped
    # gapPenalty = -5
    # similarityScore = Alignment.SimilarityScore(
    #     similarityDict=Alignment.PAM250)
    # globalAlignmentLinearGapPenalty(seq1,seq2, gapPenalty, similarityScore)

    correct = (score_res == score_sol) and is1 and is2
    return correct


def solve(seq1: str, seq2: str) -> tuple[int, str, str]:
    gapPenalty = -5
    similarityScore = Alignment.SimilarityScore(
        similarityDict=Alignment.PAM250)
    alignedSeq1, alignedSeq2, H, max_ind = Alignment.localAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)

    i, j = max_ind[0]
    max_score = H[i][j]
    return max_score, strip_inserts(alignedSeq1), strip_inserts(alignedSeq2)


def load_results(path: str) -> tuple[int, str, str]:
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
    path = f'{cwd}/rosalind_loca_1.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    gapPenalty = -5
    similarityScore = Alignment.SimilarityScore(
        similarityDict=Alignment.PAM250)
    alignedSeq1, alignedSeq2, H, max_ind = Alignment.localAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)

    i, j = max_ind[0]
    max_score = H[i][j]
    print(max_score)
    print(strip_inserts(alignedSeq1))
    print(strip_inserts(alignedSeq2))

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(max_score), 'w')
    writeTextFile(result_path, strip_inserts(alignedSeq1), 'a')
    writeTextFile(result_path, strip_inserts(alignedSeq2), 'a')

    correct = solve_and_check(path)
    print(correct)
