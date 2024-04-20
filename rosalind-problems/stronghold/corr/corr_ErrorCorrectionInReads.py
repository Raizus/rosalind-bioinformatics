from collections import defaultdict
from BioInfoToolkit.Alignment.Alignment import hamming_distance
from BioInfoToolkit.Sequences.SequenceUtils import reverseComplement

import os
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
As is the case with point mutations, the most common type of sequencing error occurs when a single nucleotide from a read is interpreted incorrectly.

    Given: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format. Some of these reads were generated with a single-nucleotide error. For each read s in the dataset, one of the following applies:

        - s was correctly sequenced and appears in the dataset at least twice (possibly as a reverse complement);

        - s is incorrect, it appears in the dataset exactly once, and its Hamming distance is 1 with respect to exactly one correct read in the dataset (or its reverse complement).

    Return: A list of all corrections in the form "[old read]->[new read]". (Each correction must be a single symbol substitution, and you may return the corrections in any order.)
"""


def match_reads(sequences: dict[str, str]):
    keys = list(sequences.keys())

    match_dict: defaultdict[str, set[str]] = defaultdict(set)
    incorrect_pairs: list[tuple[str, str]] = []
    for i, key1 in enumerate(keys[:-1]):
        seq1 = sequences[key1]
        incorrect_count = 0
        for _, key2 in enumerate(keys[i+1:], i+1):
            seq2 = sequences[key2]
            d1 = hamming_distance(seq1, seq2)
            d2 = hamming_distance(seq1, reverseComplement(seq2, 'DNA'))
            if d1 == 0 or d2 == 0:
                match_dict[key1].add(key2)
                match_dict[key2].add(key1)
            elif d1 == 1 or d2 == 1:
                incorrect_count += 1
                incorrect_pairs.append((key1, key2))

    incorrect_reads: dict[str, str] = dict()
    for (key1, key2) in incorrect_pairs:
        if key1 in match_dict and key2 not in match_dict:
            incorrect_reads[key2] = key1
        elif key2 in match_dict and key1 not in match_dict:
            incorrect_reads[key1] = key2

    return match_dict, incorrect_reads


def verify(result: list[str], solution: list[str]) -> bool:
    correct = all(r in solution for r in result) and all(r in result for r in solution)
    return correct


def solve(fasta_dict: dict[str, str]) -> list[str]:
    _, incorrect_reads = match_reads(fasta_dict)

    corrections: list[str] = []
    for incorrect_key, correct_key in incorrect_reads.items():
        if hamming_distance(fasta_dict[incorrect_key], fasta_dict[correct_key]) == 1:
            out = f"{fasta_dict[incorrect_key]}->{fasta_dict[correct_key]}"
        else:
            out = f"{fasta_dict[incorrect_key]}->{reverseComplement(fasta_dict[correct_key], 'DNA')}"
        corrections.append(out)
    return corrections


def load_results(path: str) -> list[str]:
    lines = readTextFile(path)
    corrections = [line for line in lines if len(line) and not line.isspace()]
    return corrections


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    corrections = solve(fasta_dict)

    solution = load_results(solution_path)
    correct = verify(corrections, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_corr_1.txt'

    fasta_dict = read_FASTA(path)
    corrections = solve(fasta_dict)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for correction in corrections:
        writeTextFile(result_path, correction, 'a')
        print(correction)

    correct = solve_and_check(path)
    print(correct)

