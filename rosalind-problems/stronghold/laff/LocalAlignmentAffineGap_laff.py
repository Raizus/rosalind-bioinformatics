from BioInfoToolkit import Alignment
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.SequenceUtils import is_subsequence_str
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/laff/

Given: Two protein strings s and t in FASTA format (each having length at most 10,000 aa).

Return: The maximum local alignment score of s and t, followed by substrings r and u of s and t, respectively, that correspond to the optimal local alignment of s and t. Use:
    - The BLOSUM62 scoring matrix.
    - Gap opening penalty equal to 11.
    - Gap extension penalty equal to 1.

If multiple solutions exist, then you may output any one.

"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res = result[0]
    score_sol = solution[0]

    res_seq1, res_seq2 = result[1], result[2]
    is1 = is_subsequence_str(seq1, strip_inserts(res_seq1))
    is2 = is_subsequence_str(seq2, strip_inserts(res_seq2))

    # openingPenalty: int = -11
    # extendingPenalty: int = -1
    # similarityScore = Alignment.SimilarityScore(
    #     similarityDict=Alignment.BLOSUM62)
    # actual_score = Alignment.linearGapSemiGlobalAlignmentCost(
    #     alignedSeq1, alignedSeq2, gapPenalty, similarityScore)

    correct = (score_res == score_sol) and is1 and is2
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    openingPenalty: int = -11
    extendingPenalty: int = -1
    similarityScore = Alignment.SimilarityScore(
        similarityDict=Alignment.BLOSUM62)
    alignedSeq1, alignedSeq2, V, max_ind = Alignment.smithWatermanAffineGap(
        seq1, seq2, openingPenalty, extendingPenalty, similarityScore)
    i, j = max_ind[0]
    score = V[i][j]

    return score, strip_inserts(alignedSeq1), strip_inserts(alignedSeq2)


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
    path = f'{cwd}/rosalind_laff_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1, seq2 = seqs[0], seqs[1]

    print(len(seq1))
    print(len(seq2))

    score, ss1, ss2 = solve(seq1, seq2)

    print(score)
    print(ss1)
    print(ss2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')
    writeTextFile(result_path, strip_inserts(ss1), 'a')
    writeTextFile(result_path, strip_inserts(ss2), 'a')

    # correct = solve_and_check(path)
    # print(correct)
