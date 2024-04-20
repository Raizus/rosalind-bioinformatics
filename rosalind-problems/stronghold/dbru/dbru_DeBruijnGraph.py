from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromReads
from BioInfoToolkit.Sequences.SequenceUtils import deBruijnGraph
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import re

"""
    https://rosalind.info/problems/dbru/

Consider a set S of (k+1)-mers of some unknown DNA string. Let Src denote the set containing all reverse complements of the elements of S. (recall from “Counting Subsets” that sets are not allowed to contain duplicate elements).

The de Bruijn graph Bk of order k corresponding to S∪Src is a digraph defined in the following way:

    - Nodes of Bk correspond to all k-mers that are present as a substring of a (k+1)-mer from S∪Src.
    - Edges of Bk are encoded by the (k+1)-mers of S∪Src in the following way: for each (k+1)-mer r in S∪Src, form a directed edge (r[1:k], r[2:k+1]).

    Given: A collection of up to 1000 (possibly repeating) DNA strings of equal length (not exceeding 50 bp) corresponding to a set S of (k+1)-mers.

    Return: The adjacency list corresponding to the de Bruijn graph corresponding to S∪Src.
    """


def verify(result: list[tuple[str, str]], solution: list[tuple[str, str]]) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(kmers: list[str]) -> list[tuple[str, str]]:
    k = len(kmers[0]) - 1

    graph = deBruijnGraph(kmers, k)
    edge_list: list[tuple[str,str]] = []
    for n1, n2s in graph.items():
        for n2 in n2s:
            edge_list.append((n1,n2))
    return edge_list

def parse_edge_text(text: str) -> tuple[str, str]:
    match = re.match(r'\((\w+),\s*(\w+)\)', text)
    if match is None:
        raise Exception(f"{text} does not have the correct format.")

    edge = (match[1], match[2])
    return edge


def load_results(path: str) -> list[tuple[str,str]]:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    edge_list: list[tuple[str, str]] = []
    for line in lines:
        n1, n2 = parse_edge_text(line)
        edge_list.append((n1,n2))

    return edge_list


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    kmers = [line for line in lines if len(line) and not line.isspace()]
    edge_list = solve(kmers)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(edge_list, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_dbru_1.txt'

    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]
    edge_list = solve(kmers)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for n1,n2 in edge_list:            
        out = f"({n1}, {n2})"
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
