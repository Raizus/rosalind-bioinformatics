from BioInfoToolkit import Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/sims/

Given a string s and a motif t, an alignment of a substring of s against all of t is called a fitting alignment. Our aim is to find a substring s' of s that maximizes an alignment score with respect to t.

Note that more than one such substring of s may exist, depending on the particular strings and alignment score used. One candidate for scoring function is the one derived from edit distance; In this problem, we will consider a slightly different alignment score, in which all matched symbols count as +1 and all mismatched symbols (including insertions and deletions) receive a cost of -1. Let's call this scoring function the mismatch score. See Figure 1 for a comparison of global, local, and fitting alignments with respect to mismatch score.

    Given: Two DNA strings s and t, where s has length at most 10 kbp and t represents a motif of length at most 1 kbp.

    Return: An optimal fitting alignment score with respect to the mismatch score defined above, followed by an optimal fitting alignment of a substring of s against t. If multiple such alignments exist, then you may output any one.

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
    actual_score = Alignment.linearGapGlobalAlignmentCost(
        res_seq1, res_seq2, gapPenalty, similarityScore)

    correct = (score_res == score_sol == actual_score) and is1 and is2
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty: int = -1
    similarityScore = Alignment.SimilarityScore(1, -1)
    H, max_score = Alignment.fittingAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)
    alignedSeq1, alignedSeq2 = Alignment.fittingAlignmentLinearGapBacktrack(seq1, seq2, H, gapPenalty, similarityScore)

    return max_score, alignedSeq1, alignedSeq2


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
    path = f'{cwd}/rosalind_sims_2.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1, seq2 = seqs[0], seqs[1]

    max_score, alignedSeq1, alignedSeq2 = solve(seq1, seq2)

    print(max_score)
    print(alignedSeq1)
    print(alignedSeq2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(max_score), 'w')
    writeTextFile(result_path, alignedSeq1, 'a')
    writeTextFile(result_path, alignedSeq2, 'a')

    correct = solve_and_check(path)
    print(correct)
