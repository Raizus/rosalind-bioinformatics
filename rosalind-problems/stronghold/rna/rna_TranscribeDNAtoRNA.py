from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.BioSequence import BioSequence
import os

"""
An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all occurrences of 'T' in t with 'U' in u.

    Given: A DNA string t having length at most 1000 nt.

    Return: The transcribed RNA string of t.
"""

def verify(result: str, solution: str) -> bool:
    correct = result == solution
    return correct


def solve(dna_seq: str) -> str:
    bioSeq = BioSequence(dna_seq)
    rna_seq = bioSeq.transcription()

    return rna_seq


def load_results(path: str) -> str:
    lines = readTextFile(path)
    rna_seq = lines[0]
    return rna_seq


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    result = solve(dna_seq)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_rna_1.txt'

    lines = readTextFile(path)
    dna_seq = lines[0]
    rna_seq = solve(dna_seq)
    print(rna_seq)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, rna_seq, 'w')
