from BioInfoToolkit.Sequences.SequenceUtils import reverseComplement
from BioInfoToolkit.Sequences.StringUtils import generate_d_neighborhood, hamming_distance, kmers_frequency_dictionary
from collections import Counter
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1j/

We now extend “Find the Most Frequent Words with Mismatches in a String” to find frequent words with both mismatches and reverse complements. Recall that Pattern refers to the reverse complement of Pattern.

Frequent Words with Mismatches and Reverse Complements Problem
    Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.

        Given: A DNA string Text as well as integers k and d.

        Return: All k-mers Pattern maximizing the sum Countd(Text, Pattern) + Countd(Text, Pattern) over all possible k-mers.
"""

OutputT = list[str]


def mostFrequentKmersWithMismatchesAndReverveseComplements(string: str, k: int, d: int, alphabet: set[str]):
    """Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.
    min over all k-mers of Countd(Text, Pattern) + Countd(Text, revComplementPattern)

    Args:
        string (str): _description_
        k (int): _description_
        d (int): _description_
        alphabet (set[str]): _description_

    Returns:
        _type_: _description_
    """
    kmer_freq_dict = kmers_frequency_dictionary(string, k)
    kmer_freq_with_mismatches: dict[str, int] = dict()
    for kmer, _ in kmer_freq_dict.items():
        for kmer2 in generate_d_neighborhood(kmer, d, alphabet):
            kmer2_reverse = reverseComplement(kmer2, 'DNA')
            if kmer2 in kmer_freq_with_mismatches or kmer2_reverse in kmer_freq_with_mismatches:
                continue

            freq = 0
            for kmer3, count2 in kmer_freq_dict.items():
                dist1 = hamming_distance(kmer2, kmer3)
                if dist1 <= d:
                    freq += count2
                dist2 = hamming_distance(kmer2_reverse, kmer3)
                if dist2 <= d:
                    freq += count2

            kmer_freq_with_mismatches[kmer2] = freq
            kmer_freq_with_mismatches[kmer2_reverse] = freq

    max_freq = max(kmer_freq_with_mismatches.values())
    kmers = [kmer for kmer, freq in kmer_freq_with_mismatches.items()
             if freq == max_freq]
    return kmers


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(sequence: str, k: int, d: int) -> OutputT:
    alphabet = {'A', 'C', 'G', 'T'}
    kmers = mostFrequentKmersWithMismatchesAndReverveseComplements(
        sequence, k, d, alphabet)
    return kmers


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    kmers = lines[0].split()
    return kmers


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    sequence = lines[0]
    k, d = [int(v) for v in lines[1].split()]

    result = solve(sequence, k, d)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1j_1.txt'

    lines = readTextFile(path)
    sequence = lines[0]
    k, d = [int(v) for v in lines[1].split()]

    alphabet = {'A', 'C', 'G', 'T'}
    kmers = solve(sequence, k, d)

    out = ' '.join(kmers)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

