from BioInfoToolkit.Sequences.SequenceUtils import failure_array
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/kmp/

A prefix of a length n string s is a substring s[1:j]; a suffix of s is a substring s[k:n].

The failure array of s is an array P of length n for which P[k] is the length of the longest substring s[j:k] that is equal to some prefix s[1:k-j+1], where j cannot equal 1 (otherwise, P[k] would always equal k). By convention, P[1]=0.

    Given: A DNA string s (of length at most 100 kbp) in FASTA format.

    Return: The failure array of s.
"""


def verify(result: list[int], solution: list[int]) -> bool:
    correct = len(result) == len(solution) and all(
        v1 == v2 for v1, v2 in zip(result, solution))
    return correct


def solve(seq: str) -> list[int]:
    failure = failure_array(seq)
    return failure


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    failure = [int(v) for v in lines[0].split()]
    return failure


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]
    result = solve(seq)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_kmp_1.txt'

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]
    failure = failure_array(seq)

    out = " ".join(str(val) for val in failure)

    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
