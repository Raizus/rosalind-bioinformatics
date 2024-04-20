from collections import Counter
from BioInfoToolkit.Sequences.StringUtils import mostFrequentKmersWithMismatches
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1i/

We defined a mismatch in “Compute the Hamming Distance Between Two Strings”. We now generalize “Find the Most Frequent Words in a String” to incorporate mismatches as well.

Given strings Text and Pattern as well as an integer d, we define Countd(Text, Pattern) as the total number of occurrences of Pattern in Text with at most d mismatches. For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA, AAACA, and AAAGA. Note that two of these occurrences overlap.

A most frequent k-mer with up to d mismatches in Text is simply a string Pattern maximizing Countd(Text, Pattern) among all k-mers. Note that Pattern does not need to actually appear as a substring of Text; for example, AAAAA is the most frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG, even though AAAAA does not appear exactly in this string. Keep this in mind while solving the following problem.

Frequent Words with Mismatches Problem
    Find the most frequent k-mers with mismatches in a string.

        Given: A string Text as well as integers k and d.

        Return: All most frequent k-mers with up to d mismatches in Text.
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(sequence: str, k: int, d: int) -> OutputT:
    alphabet = {'A', 'C', 'G', 'T'}
    kmers = mostFrequentKmersWithMismatches(sequence, k, d, alphabet)
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
    path = f'{cwd}/rosalind_ba1i_2.txt'

    lines = readTextFile(path)
    sequence = lines[0]
    k, d = [int(v) for v in lines[1].split()]
    kmers = solve(sequence, k, d)

    out = ' '.join(kmers)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
