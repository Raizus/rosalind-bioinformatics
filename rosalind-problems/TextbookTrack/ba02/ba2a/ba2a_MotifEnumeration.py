from BioInfoToolkit.Sequences.StringUtils import motifEnumeration
from collections import Counter
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba2a/

Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif if it appears in every string from Dna with at most d mismatches. The following algorithm finds (k,d)-motifs.

    MOTIFENUMERATION(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern in Dna
            for each k-mer Pattern' differing from Pattern by at most d
              mismatches
                if Pattern' appears in each string from Dna with at most d
                mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns

Implanted Motif Problem
    Implement MotifEnumeration (shown above) to find all (k, d)-motifs in a collection of strings.

        Given: Integers k and d, followed by a collection of strings Dna.

        Return: All (k, d)-motifs in Dna.
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(sequences: list[str], k: int, d: int) -> OutputT:
    alphabet = set(['A', 'C', 'G', 'T'])
    patterns = motifEnumeration(sequences, k, d, alphabet)
    patterns = sorted(list(patterns))
    return patterns


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    motifs = lines[0].split()
    return motifs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, d = [int(a) for a in lines[0].split()]
    sequences = lines[1:]

    result = solve(sequences, k, d)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba2a_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, d = [int(a) for a in lines[0].split()]
    sequences = lines[1:]

    strings = solve(sequences, k, d)

    out = ' '.join(strings)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
