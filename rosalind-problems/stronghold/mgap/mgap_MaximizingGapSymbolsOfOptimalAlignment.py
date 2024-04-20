from BioInfoToolkit.Sequences.SequenceUtils import longestCommonSubsequenceTable
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/mgap/

For the computation of an alignment score generalizing the edit alignment score, let m denote the score assigned to matched symbols, d denote the score assigned to mismatched non-gap symbols, and g denote the score assigned a symbol matched to a gap symbol '-' (i.e., g is a linear gap penalty).

    Given: Two DNA strings s and t in FASTA format (each of length at most 5000 bp).

    Return: The maximum number of gap symbols that can appear in any maximum score alignment of s and t with score parameters satisfying m>0, d<0, and g<0.

"""


def maximumGaps(s: str, t: str) -> int:
    C = longestCommonSubsequenceTable(s, t)
    val = C[-1][-1]

    count = len(s) + len(t) - 2*val
    return count


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq1: str, seq2: str) -> int:
    count = maximumGaps(seq1, seq2)
    return count


def load_results(path: str) -> int:
    lines = readTextFile(path)
    gap_num = int(lines[0])

    return gap_num


def solve_and_check(input_path: str) -> bool:
    FASTAdict = read_FASTA(input_path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    result = solve(seq1, seq2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_mgap_1.txt'

    FASTAdict = read_FASTA(path)
    seqs = list(FASTAdict.values())
    seq1, seq2 = seqs[0], seqs[1]

    count = maximumGaps(seq1, seq2)
    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    # correct = solve_and_check(path)
    # print(correct)

