from BioInfoToolkit import Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/smgb/

A semiglobal alignment of strings s and t is an alignment in which any gaps appearing as prefixes or suffixes of s and t do not contribute to the alignment score.

Semiglobal alignment has sometimes also been called "overlap alignment". Rosalind defines overlap alignment differently (see “Overlap Alignment”). 

    Given: Two DNA strings s and t in FASTA format, each having length at most 10 kbp.

    Return: The maximum semiglobal alignment score of s and t, followed by an alignment of s and t achieving this maximum score. Use an alignment score in which matching symbols count +1, substitutions count -1, and there is a linear gap penalty of 1. If multiple optimal alignments exist, then you may return any one.

"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res = result[0]
    score_sol = solution[0]

    res_seq1, res_seq2 = result[1], result[2]
    is1 = is_subsequence_str(seq1, strip_inserts(res_seq1))
    is2 = is_subsequence_str(seq2, strip_inserts(res_seq2))

    gapPenalty: int = -1
    similarityScore = Alignment.SimilarityScore(1, -1)
    actual_score = Alignment.linearGapSemiGlobalAlignmentCost(
        alignedSeq1, alignedSeq2, gapPenalty, similarityScore)

    correct = (score_res == score_sol == actual_score) and is1 and is2
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty: int = -1
    similarityScore = Alignment.SimilarityScore(1, -1)
    alignedSeq1, alignedSeq2, _, score = Alignment.semiglobalAligmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)

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
    path = f'{cwd}/rosalind_smgb_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1, seq2 = seqs[0], seqs[1]

    print(len(seq1))
    print(len(seq2))

    gapPenalty: int = -1
    similarityScore = Alignment.SimilarityScore(1, -1)
    alignedSeq1, alignedSeq2, _, score = Alignment.semiglobalAligmentLinearGap(seq1, seq2, gapPenalty, similarityScore)

    print(score)
    print(alignedSeq1)
    print(alignedSeq2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')
    writeTextFile(result_path, alignedSeq1, 'a')
    writeTextFile(result_path, alignedSeq2, 'a')

    # correct = solve_and_check(path)
    # print(correct)

