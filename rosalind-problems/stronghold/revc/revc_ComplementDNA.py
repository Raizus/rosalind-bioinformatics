from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.BioSequence import BioSequence
import os

"""
In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The reverse complement of a DNA string s is the string s^c formed by reversing the symbols of s, then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

    Given: A DNA string s of length at most 1000 bp.
    
    Return: The reverse complement s^c of s.
"""


def verify(result: str, solution: str) -> bool:
    correct = result == solution
    return correct


def solve(dna_seq: str) -> str:
    bioSeq = BioSequence(dna_seq)
    rev_complement = bioSeq.reverseComplement()

    return rev_complement


def load_results(path: str) -> str:
    lines = readTextFile(path)
    seq = lines[0]
    return seq


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    result = solve(dna_seq)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_revc_1.txt'

    lines = readTextFile(path)
    dna_seq = lines[0]
    rev_complement = solve(dna_seq)
    print(rev_complement)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, rev_complement, 'w')
