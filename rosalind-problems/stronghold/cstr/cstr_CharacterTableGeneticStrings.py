from collections import Counter
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/cstr/

A collection of strings is characterizable if there are at most two possible choices for the symbol at each position of the strings.

    Given: A collection of at most 100 characterizable DNA strings, each of length at most 300 bp.

    Return: A character table for which each nontrivial character encodes the symbol choice at a single position of the strings. (Note: the choice of assigning '1' and '0' to the two states of each SNP in the strings is arbitrary.)
"""


def verify(result: list[str], solution: list[str]) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(taxa: list[str]) -> list[str]:
    result: list[str] = []
    for i, ns in enumerate(zip(*taxa)):
        counts = Counter(ns)
        # the split is characterizable if it has at most two possible nucleotides at each position
        nn = len(counts)
        if nn < 2:
            continue
        # this is a split
        elif nn == 2 and all(count > 1 for count in counts.values()):
            character_table = [int(taxon[i] == ns[0]) for taxon in taxa]
            result.append(''.join(str(c) for c in character_table))
    return result


def load_results(path: str) -> list[str]:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    return lines


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    taxa = [line for line in lines if len(line) and not line.isspace()]

    character_table = solve(taxa)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(character_table, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_cstr_1.txt'

    lines = readTextFile(path)
    taxa = [line for line in lines if len(line) and not line.isspace()]

    character_table = solve(taxa)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for c in character_table:
        print(c)
        writeTextFile(result_path, c, 'a')

    correct = solve_and_check(path)
    print(correct)
