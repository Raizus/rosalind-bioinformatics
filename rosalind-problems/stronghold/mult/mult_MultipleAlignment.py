from BioInfoToolkit.Alignment.MultiAlignment import getMultiAlignScoringFunc, multiAlignmentCost, multisequenceAlign
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Sequences.StringUtils import strip_inserts

"""
https://rosalind.info/problems/mult/

A multiple alignment of a collection of three or more strings is formed by adding gap symbols to the strings to produce a collection of augmented strings all having the same length.

A multiple alignment score is obtained by taking the sum of an alignment score over all possible pairs of augmented strings. The only difference in scoring the alignment of two strings is that two gap symbols may be aligned for a given pair (requiring us to specify a score for matched gap symbols).

    Given: A collection of four DNA strings of length at most 10 bp in FASTA format.

    Return: A multiple alignment of the strings having maximum score, where we score matched symbols 0 (including matched gap symbols) and all mismatched symbols -1 (thus incorporating a linear gap penalty of 1).

"""


def verify(result: tuple[int, list[str]], solution: tuple[int, list[str]], seqs: list[str]) -> bool:
    score_res = result[0]
    score_sol = solution[0]

    aligned_seqs_res = result[1]

    scoring_func = getMultiAlignScoringFunc(0, -1)
    alignment_score = multiAlignmentCost(aligned_seqs_res, scoring_func)
    are_subseqs = all(seq == strip_inserts(aligned_seq) for seq, aligned_seq in
                      zip(seqs, aligned_seqs_res))

    correct = (score_res == score_sol == alignment_score) and are_subseqs
    return correct


def solve(seqs: list[str]) -> tuple[int, list[str]]:
    scoringFunc = getMultiAlignScoringFunc(0, -1)
    aligned_seqs, final_score = multisequenceAlign(seqs, scoringFunc)
    return final_score, aligned_seqs


def load_results(path: str) -> tuple[int, list[str]]:
    lines = readTextFile(path)
    score = int(lines[0])
    aligned_seqs = lines[1:]
    aligned_seqs = [s for s in aligned_seqs if len(s) and not s.isspace()]

    return score, aligned_seqs


def solve_and_check(input_path: str) -> bool:
    fasta_dict = read_FASTA(input_path)
    seqs = list(fasta_dict.values())

    result = solve(seqs)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, seqs)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_mult_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())

    scoring_func = getMultiAlignScoringFunc(0, -1)
    aligned_seqs, final_score = multisequenceAlign(
        seqs, scoring_func)

    print(final_score)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(final_score), 'w')
    for seq in aligned_seqs:
        print(seq)
        writeTextFile(result_path, seq, 'a')

    # correct = solve_and_check(path)
    # print(correct)
