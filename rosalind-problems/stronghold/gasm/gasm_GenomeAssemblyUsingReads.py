from collections import Counter, defaultdict
from BioInfoToolkit.IO.IO import readTextFile
from BioInfoToolkit.Sequences.SequenceUtils import deBruijnGraph, reverseComplement
from BioInfoToolkit.Sequences.StringUtils import kmer_gen, rotationally_equivalent

from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromReads
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
    https://rosalind.info/problems/pcov/

A circular string is a string that does not have an initial or terminal element; instead, the string is viewed as a necklace of symbols. We can represent a circular string as a string enclosed in parentheses. For example, consider the circular DNA string (ACGTAC), and note that because the string "wraps around" at the end, this circular string can equally be represented by (CGTACA), (GTACAC), (TACACG), (ACACGT), and (CACGTA). The definitions of substrings and superstrings are easy to generalize to the case of circular strings (keeping in mind that substrings are allowed to wrap around).

    Given: A collection of (error-free) DNA k-mers (k≤50) taken from the same strand of a circular chromosome. In this dataset, all k-mers from this strand of the chromosome are present, and their de Bruijn graph consists of exactly one simple cycle.

    Return: A cyclic superstring of minimal length containing the reads (thus corresponding to a candidate cyclic chromosome).
"""

def generate_all_kmers(seqs: list[str], k: int):
    for seq in seqs:
        for kmer in kmer_gen(seq, k):
            yield kmer


def find_cycle(graph: defaultdict[str, set[str]]):
    startNode = next(iter(graph))
    n1 = startNode
    n2: str | None = None

    path: list[str] = []
    path.append(n1)
    while n2 != startNode:
        n2s = graph.get(n1)
        if not n2s:
            return []
        n2 = next(iter(n2s))
        if (n2 != startNode):
            path.append(n2)

        n1 = n2
    return path


def verify(result: str, solution: str, reads: list[str]) -> bool:
    # TODO: maybe they're not necesarily equivalent and reads should be used to check the result
    rot_equiv = rotationally_equivalent(result, solution) or \
        rotationally_equivalent(result, reverseComplement(solution, 'DNA'))
    correct = len(result) == len(solution) and rot_equiv
    return correct


def solve(kmers: list[str]) -> str:
    max_k = len(kmers[0])
    superstring = ""

    for k1 in range(max_k, 2, -1):
        kmers = list(generate_all_kmers(kmers, k1))
        graph = deBruijnGraph(kmers, k1-1)

        # check if it has exactly 2 cycles
        path = find_cycle(graph)

        # if so reconstruct the superstring
        if len(path):
            superstring = path[0]
            for seq in path[1:-k1+2]:
                superstring += seq[-1]
            break
    return superstring


def load_results(path: str) -> str:
    lines = readTextFile(path)
    string = lines[0]
    return string


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    kmers = [line for line in lines if len(line) and not line.isspace()]
    string = solve(kmers)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(string, solution, kmers)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_gasm_1.txt'

    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    superstring = solve(kmers)

    print(superstring)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, superstring, 'w')

    correct = solve_and_check(path)
    print(correct)
