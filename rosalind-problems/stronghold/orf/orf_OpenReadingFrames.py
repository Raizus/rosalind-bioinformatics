from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.BioSequence import BioSequence
import os

"""
Either strand of a DNA double helix can serve as the coding strand for RNA transcription. Hence, a given DNA string implies six total reading frames, or ways in which the same region of DNA can be translated into amino acids: three reading frames result from reading the string itself, whereas three more result from reading its reverse complement.

An open reading frame (ORF) is one which starts from the start codon and ends by stop codon, without any other stop codons in between. Thus, a candidate protein string is derived by translating an open reading frame into amino acids until a stop codon is reached.

    Given: A DNA string s of length at most 1 kbp in FASTA format.

    Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.
"""


def verify(result: list[str], solution: list[str]) -> bool:
    return set(result) == set(solution)


def solve(seq: str) -> list[str]:
    bio_seq = BioSequence(seq)
    proteins = bio_seq.all_proteins_from_rfs(ordered=True)
    return proteins


def load_results(path: str) -> list[str]:
    proteins = readTextFile(path)
    proteins = [prot for prot in proteins if not prot.isspace() and len(prot)]
    return proteins


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
    path = f'{cwd}/rosalind_orf_1.txt'

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]
    bio_seq = BioSequence(seq)
    proteins = bio_seq.all_proteins_from_rfs(ordered=True)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, '\n'.join(proteins), 'w')
    for p in proteins:
        print(p)
