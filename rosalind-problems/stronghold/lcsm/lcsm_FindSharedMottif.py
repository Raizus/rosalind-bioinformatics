from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

from BioInfoToolkit.Sequences.SequenceUtils import longestCommonString

"""
A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".

Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

    Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.

    Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)
"""


def verify(fasta_dict: dict[str, str], result: str, solution: str) -> bool:
    len_equals = len(result) == len(solution)
    seqs = list(fasta_dict.values())

    correct = len_equals and all(solution in seq for seq in seqs)
    return correct


def solve(fasta_dict: dict[str, str]) -> str:
    seqs = list(fasta_dict.values())
    lcs = longestCommonString(seqs)
    return lcs


def load_results(path: str) -> str:
    lines = readTextFile(path)
    lcs = lines[0]
    return lcs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(input_path)
    result = solve(fasta_dict)
    solution = load_results(solution_path)

    correct = verify(fasta_dict, result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_lcsm_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    lcs = longestCommonString(seqs)

    print(lcs)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, lcs, 'w')

