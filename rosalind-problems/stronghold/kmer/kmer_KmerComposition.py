from itertools import product
from BioInfoToolkit.Sequences.StringUtils import kmers_frequency_dictionary
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.

    Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.

    Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)
"""


def verify(result: list[int], solution: list[int]) -> bool:
    correct = len(result) == len(solution) and all(v1 == v2 for v1, v2 in zip(result, solution))
    return correct


def solve(seq: str, kmer_l: int) -> list[int]:
    kmer_freq = kmers_frequency_dictionary(seq, kmer_len)
    # add remaining kmers (those not appearing in the sequence)
    for p in product("ACTG", repeat=kmer_len):
        kmer_freq.setdefault("".join(str(n) for n in p), 0)

    sorted_keys = sorted(kmer_freq.keys())
    freqs = [kmer_freq[key] for key in sorted_keys]

    return freqs


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    freqs = [int(v) for v in lines[0].split()]
    return freqs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    kmer_len = 4
    seq = list(fasta_dict.values())[0]
    result = solve(seq, kmer_len)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_kmer_1.txt'

    fasta_dict = read_FASTA(path)
    kmer_len = 4
    seq = list(fasta_dict.values())[0]
    freqs = solve(seq, kmer_len)

    out = " ".join(str(f) for f in freqs)

    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    # correct = solve_and_check(path)
    # print(correct)
