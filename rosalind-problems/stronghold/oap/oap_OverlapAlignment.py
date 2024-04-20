import BioInfoToolkit.Alignment as Alignment
from BioInfoToolkit.Alignment.Alignment import linearGapOverlapAlignmentCost
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/oap/

An overlap alignment between two strings s and t is a local alignment of a suffix of s with a prefix of t. An optimal overlap alignment will therefore maximize an alignment score over all such substrings of s and t.

The term "overlap alignment" has also been used to describe what Rosalind defines as a semiglobal alignment. See “Semiglobal Alignment” for details.

    Given: Two DNA strings s and t in FASTA format, each having length at most 10 kbp.

    Return: The score of an optimal overlap alignment of s
and t, followed by an alignment of a suffix s' of s and a prefix t' of t achieving this optimal score. Use an alignment score in which matching symbols count +1, substitutions count -2, and there is a linear gap penalty of 2. If multiple optimal alignments exist, then you may return any one.

"""

OutputT = tuple[int, str, str]


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    score_res = result[0]
    score_sol = solution[0]

    res_seq1, res_seq2 = result[1], result[2]
    is1 = strip_inserts(res_seq1) == seq1[-len(strip_inserts(res_seq1)):]
    is2 = strip_inserts(res_seq2) == seq2[:len(strip_inserts(res_seq2))]

    gapPenalty: int = -2
    similarityScore = Alignment.SimilarityScore(1, -2)
    actual_score = linearGapOverlapAlignmentCost(res_seq1, res_seq2, gapPenalty, similarityScore)

    correct = (score_res == score_sol == actual_score) and is1 and is2
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty: int = -2
    similarityScore = Alignment.SimilarityScore(1, -2)
    aligned_s, aligned_t, _, score = Alignment.overlapAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)

    new_s, new_t = trim_s_t(aligned_s, aligned_t)
    return score, new_s, new_t


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


def trim_s_t(s: str, t: str) -> tuple[str, str]:
    m = len(s)
    new_s = s.rstrip('-')
    ds = m - len(new_s)
    t = t[:-ds]
    n = len(t)

    new_t = t.lstrip('-')
    dt = n - len(new_t)
    new_s = new_s[dt:]

    return new_s, new_t


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_oap_0.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    # print(len(seq1))
    # print(len(seq2))

    gapPenalty: int = -2
    similarityScore = Alignment.SimilarityScore(1, -2)
    aligned_s, aligned_t, _, score = Alignment.overlapAlignmentLinearGap(
        seq1, seq2, gapPenalty, similarityScore)
    
    score, aligned_s, aligned_t = solve(seq1, seq2)

    print(score)
    print(aligned_s)
    print(aligned_t)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')
    writeTextFile(result_path, aligned_s, 'a')
    writeTextFile(result_path, aligned_t, 'a')

    correct = solve_and_check(path)
    print(correct)

