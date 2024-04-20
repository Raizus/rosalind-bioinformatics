from BioInfoToolkit.Phylogeny.Phylo import find_quartets
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import re

"""
https://rosalind.info/problems/qrt/

A partial split of a set S of n taxa models a partial character and is denoted by A|B, where A and B are still the two disjoint subsets of taxa divided by the character. Unlike in the case of splits, we do not necessarily require that A∪B=S; (A∪B)c corresponds to those taxa for which we lack conclusive evidence regarding the character.

We can assemble a collection of partial characters into a generalized partial character table C
in which the symbol x is placed in Ci,j if we do not have conclusive evidence regarding the jth taxon with respect to the ith partial character.

A quartet is a partial split A|B in which both A and B contain precisely two elements. For the sake of simplicity, we often will consider quartets instead of partial characters. We say that a quartet A|B is inferred from a partial split C|D if A⊆C and B⊆D (or equivalently A⊆D and B⊆C). For example, {1,3}|{2,4} and {3,5}|{2,4} can be inferred from {1,3,5}|{2,4}.

    Given: A partial character table C.

    Return: The collection of all quartets that can be inferred from the splits corresponding to the underlying characters of C.
"""


def format_pair(pair: tuple[str, str]) -> str:
    res = "{" + f"{pair[0]}, {pair[1]}" + "}"
    return res


def pair_idx_to_pair(taxa: list[str], pair: tuple[int, int]) -> tuple[str, str]:
    return (taxa[pair[0]], taxa[pair[1]])


def quartets_are_equal(q1: tuple[tuple[str, str], tuple[str, str]], 
                       q2: tuple[tuple[str, str], tuple[str, str]]):
    p1, p2 = q1
    p3, p4 = q2

    p1r = (p1[1], p1[0])
    p2r = (p2[1], p2[0])
    
    p1_in_q2 = p1 == p3 or p1 == p4 or p1r == p3 or p1r == p4
    p2_in_q2 = p2 == p3 or p2 == p4 or p2r == p3 or p2r == p4

    return p1_in_q2 and p2_in_q2


def verify(result: list[tuple[tuple[str, str], tuple[str, str]]], 
           solution: list[tuple[tuple[str, str], tuple[str, str]]]) -> bool:
    if len(result) != len(solution):
        return False
    for q1 in result:
        if not any(quartets_are_equal(q1, q2) for q2 in solution):
            return False
    return True


def solve(character_table: list[str], taxa: list[str]) -> list[tuple[tuple[str, str], tuple[str, str]]]:
    quartets_i = find_quartets(character_table, True)

    quartets: list[tuple[tuple[str, str], tuple[str, str]]] = []
    for p1_i, p2_i in quartets_i:
        p1 = pair_idx_to_pair(taxa, p1_i)
        p2 = pair_idx_to_pair(taxa, p2_i)
        q = (p1, p2)
        quartets.append(q)

    return quartets


def parse_quartets(string: str) -> tuple[tuple[str, str], tuple[str, str]]:
    match = re.match(r'\{(\w+),\s*(\w+)\}\s*\{(\w+),\s*(\w+)\}', string)
    if match is None:
        raise Exception(f"{string} does not have the correct quartet format.")
    quartet = ((match[1], match[2]), (match[3], match[4]))

    return quartet


def load_results(path: str) -> list[tuple[tuple[str, str], tuple[str, str]]]:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    quartets: list[tuple[tuple[str,str],tuple[str,str]]] = []
    for line in lines:
        quartet = parse_quartets(line)
        quartets.append(quartet)
    return quartets


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    taxa = lines[0].split(' ')
    character_table = lines[1:]

    quartets = solve(character_table, taxa)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(quartets, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_qrt_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    taxa = lines[0].split(' ')
    character_table = lines[1:]

    quartets = solve(character_table, taxa)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for p0, p1 in quartets:
        text = f"{format_pair(p0)} {format_pair(p1)}"
        writeTextFile(result_path, text, 'a')
        print(text)

    correct = solve_and_check(path)
    print(correct)
