from BioInfoToolkit.Sequences.BioSequence import BioSequence
from BioInfoToolkit.IO import read_FASTA

from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.

    Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.

    Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)
"""


def verify(result: str, solution: str) -> bool:
    correct = result == solution
    return correct


def solve(seq: str, introns: list[str]) -> str:
    bio_seq = BioSequence(seq)
    prot = bio_seq.RNASplicing(introns)
    return prot


def load_results(path: str) -> str:
    lines = readTextFile(path)
    prot = lines[0]
    return prot


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seqs = sorted(fasta_dict.values(), key=lambda s: len(s), reverse=True)
    seq = seqs[0]
    introns = seqs[1:]

    prot = solve(seq, introns)

    solution = load_results(solution_path)

    correct = verify(prot, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_splc_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = sorted(fasta_dict.values(), key=lambda s: len(s), reverse=True)
    seq = seqs[0]
    introns = seqs[1:]

    bio_seq = BioSequence(seq)
    prot = bio_seq.RNASplicing(introns)

    print(prot)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, prot, 'w')

    correct = solve_and_check(path)
    print(correct)
