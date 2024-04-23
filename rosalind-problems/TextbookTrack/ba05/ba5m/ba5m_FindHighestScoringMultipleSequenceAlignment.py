from BioInfoToolkit.Alignment.MultiAlignment import multiAlignmentCost, multisequenceAlign
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/ba5m/

Multiple Longest Common Subsequence Problem
    Find a longest common subsequence of multiple strings.

        Given: Three DNA strings.

        Return: The maximum score of a multiple alignment of these three strings, followed by a multiple alignment of the three strings achieving this maximum. Use a scoring function in which the score of an alignment column is 1 if all three symbols are identical and 0 otherwise. (If more than one multiple alignment achieve the maximum, you may return any one.)
"""

OutputT = tuple[int, list[str]]


def scoringFunc(chars: list[str]) -> int:
    all_equal = all(chars[0] == char for char in chars[1:])
    totalScore = 1 if all_equal else 0
    return totalScore


def verify(result: OutputT, solution: OutputT, seqs: list[str]) -> bool:
    score_res = result[0]
    score_sol = solution[0]

    aligned_seqs_res = result[1]

    alignment_score = multiAlignmentCost(aligned_seqs_res, scoringFunc)
    are_subseqs = all(seq == strip_inserts(aligned_seq) for seq, aligned_seq in
                      zip(seqs, aligned_seqs_res))

    correct = (score_res == score_sol == alignment_score) and are_subseqs

    return correct


def solve(seqs: list[str]) -> OutputT:
    aligned_seqs, final_score = multisequenceAlign(
        seqs, scoringFunc)

    return final_score, aligned_seqs


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    score = int(lines[0])
    seqs = [line for line in lines[1:] if len(line) and not line.isspace()]

    return score, seqs


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seqs = [line for line in lines if len(line) and not line.isspace()]

    result = solve(seqs)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, seqs)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba5m_1.txt'

    lines = readTextFile(path)
    seqs = [line for line in lines if len(line) and not line.isspace()]

    score, aligned_seqs = solve(seqs)

    print(score)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')
    for seq in aligned_seqs:
        print(seq)
        writeTextFile(result_path, seq, 'a')

    correct = solve_and_check(path)
    print(correct)

