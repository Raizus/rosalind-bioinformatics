from BioInfoToolkit import Alignment
from BioInfoToolkit.Alignment.Alignment import SimilarityScore
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/osym/

Say that we have two strings s and t of respective lengths m and n and an alignment score. Let's define a matrix M corresponding to s and t by setting Mj,k equal to the maximum score of any alignment that aligns s[j] with t[k]. So each entry in M can be equal to at most the maximum score of any alignment of s and t.

    Given: Two DNA strings s and t in FASTA format, each having length at most 1000 bp.

    Return: The maximum alignment score of a global alignment of s and t, followed by the sum of all elements of the matrix M corresponding to s and t that was defined above. Apply the mismatch score introduced in “Finding a Motif with Modifications”.

"""

OutputT = tuple[int, int]


def matchSymbolCostMatrix(seq1: str, seq2: str, gap_penalty: int, similarity_score: SimilarityScore):
    m = len(seq1)
    n = len(seq2)

    left_cost_gen = Alignment.globalAlignmentLinearGapGen(
        seq1, seq2, gap_penalty, similarity_score)
    F2, _ = Alignment.globalAlignmentLinearGapPenalty(
        seq1[::-1], seq2[::-1], gap_penalty, similarity_score)

    M: list[list[int]] = [[0 for _ in range(n)] for _ in range(m)]
    for i, (c1, left_cost_line) in enumerate(zip(seq1, left_cost_gen)):
        for j, c2 in enumerate(seq2):
            left_cost = left_cost_line[j]
            right_cost = F2[m-i-1][n-j-1]
            forced_match_cost = similarity_score.score(c1, c2)
            total_cost = left_cost + forced_match_cost + right_cost
            M[i][j] = total_cost

    return M


def verify(result: OutputT, solution: OutputT, seq1: str, seq2: str) -> bool:
    correct = (result == solution)
    return correct


def solve(seq1: str, seq2: str) -> OutputT:
    gapPenalty = -1
    matchScore = 1
    mismatchScore = -1
    similarityScore = Alignment.SimilarityScore(matchScore, mismatchScore)

    M = matchSymbolCostMatrix(seq1, seq2, gapPenalty, similarityScore)

    total_sum = sum([sum(l) for l in M])
    max_val = max([max(l) for l in M])

    return max_val, total_sum


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    score = int(lines[0])
    sum_val = int(lines[1])

    return score, sum_val


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
    path = f'{cwd}/rosalind_osym_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1, seq2 = seqs[0], seqs[1]

    # print(len(seq1))
    # print(len(seq2))

    max_val, sum_val = solve(seq1, seq2)

    print(max_val)
    print(sum_val)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(max_val), 'w')
    writeTextFile(result_path, str(sum_val), 'a')

    correct = solve_and_check(path)
    print(correct)

