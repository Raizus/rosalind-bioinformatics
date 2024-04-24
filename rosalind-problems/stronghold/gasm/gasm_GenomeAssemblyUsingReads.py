from collections import defaultdict
from BioInfoToolkit.Sequences.SequenceUtils import deBruijnGraph, reverseComplement
from BioInfoToolkit.Sequences.StringUtils import kmer_gen, rotationally_equivalent
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/gasm/

Problem

A directed cycle is simply a cycle in a directed graph in which the head of one edge is equal to the tail of the next (so that every edge in the cycle is traversed in the same direction).

For a set of DNA strings S and a positive integer k, let Sk denote the collection of all possible k-mers of the strings in S.

    Given: A collection S of (error-free) reads of equal length (not exceeding 50 bp). In this dataset, for some positive integer k, the de Bruijn graph Bk on Sk+1 ∪ Srck+1 consists of exactly two directed cycles.

    Return: A cyclic superstring of minimal length containing every read or its reverse complement.
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
    rot_equiv = rotationally_equivalent(result, solution) or \
        rotationally_equivalent(result, reverseComplement(solution, 'DNA'))
    
    # check if the superstring contains all the reads or their reverse complement
    k = len(reads[0])
    idxs: set[int] = set()
    for kmer in kmer_gen(result, k, True):
        idxs_ = [i for i, read in enumerate(
            reads) if kmer == read or kmer == reverseComplement(read, 'DNA')]
        idxs.update(idxs_)

    correct = len(result) == len(solution) and rot_equiv and len(idxs) == len(reads)
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
