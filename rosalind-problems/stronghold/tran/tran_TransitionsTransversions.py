from BioInfoToolkit.Sequences.SequenceUtils import count_transitions_transversions
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
For DNA strings s1 and s2 having the same length, their transition/transversion ratio R(s1,s2) is the ratio of the total number of transitions to the total number of transversions, where symbol substitutions are inferred from mismatched corresponding symbols as when calculating Hamming distance (see “Counting Point Mutations”).

    Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

    Return: The transition/transversion ratio R(s1,s2).
"""


def verify(result: float, solution: float) -> bool:
    error = abs(result-solution)
    correct = error <= 0.001
    return correct


def solve(seq1: str, seq2: str) -> float:
    (transitions, tranversions) = count_transitions_transversions(seq1, seq2)
    ratio = transitions/tranversions
    return ratio


def load_results(path: str) -> float:
    lines = readTextFile(path)
    ratio = float(lines[0])
    return ratio


def solve_and_check(input_path: str) -> bool:
    fasta_dict = read_FASTA(input_path)
    seqs = list(fasta_dict.values())
    seq1 = seqs[0]
    seq2 = seqs[1]

    result = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_tran_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    seq1 = seqs[0]
    seq2 = seqs[1]

    ratio = solve(seq1, seq2)

    print(ratio)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(ratio), 'w')

    correct = solve_and_check(path)
    print(correct)
