from collections import Counter
from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromReads
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import re

"""
https://rosalind.info/problems/ba3e/

Given an arbitrary collection of k-mers Patterns (where some k-mers may appear multiple times), we define CompositionGraph(Patterns) as a graph with |Patterns| isolated edges. Every edge is labeled by a k-mer from Patterns, and the starting and ending nodes of an edge are labeled by the prefix and suffix of the k-mer labeling that edge. We then define the de Bruijn graph of Patterns, denoted DeBruijn(Patterns), by gluing identically labeled nodes in CompositionGraph(Patterns), which yields the following algorithm.

    DEBRUIJN(Patterns)
        represent every k-mer in Patterns as an isolated edge between its prefix and suffix
        glue all nodes with identical labels, yielding the graph DeBruijn(Patterns)
        return DeBruijn(Patterns) 

De Bruijn Graph from k-mers Problem
    Construct the de Bruijn graph from a collection of k-mers.

        Given: A collection of k-mers Patterns.

        Return: The de Bruijn graph DeBruijn(Patterns), in the form of an adjacency list.        
"""

OutputT = dict[str, list[str]]


def verify(result: OutputT, solution: OutputT) -> bool:
    if set(result.keys()) != set(solution.keys()):
        return False

    for key1, val1 in result.items():
        if Counter(solution[key1]) != Counter(val1):
            return False
    return True


def solve(kmers: list[str]) -> OutputT:
    k = len(kmers[0])
    _g = deBruijnMultiGraphFromReads(kmers, k)
    g = DeBruijnMultiGraph(_g, k)
    # dot = g.draw_dot()
    # dot.view()

    # string = g.reconstructStringFromEulerianPath()
    # print(string)

    sorted_nodes = sorted([nodeId for nodeId in g.graph.nodes])
    edge_dict: OutputT = dict()
    for node1 in sorted_nodes:
        adjacency_dict = g.graph.adj[node1]
        adj_seqs: list[str] = []
        for node2, edges in adjacency_dict.items():
            adj_seqs.extend([node2]*len(edges))

        if len(adj_seqs):
            edge_dict[node1] = adj_seqs

    return edge_dict


def parse_line(line: str) -> tuple[str, list[str]]:
    match = re.match(r'(\w+)\s?->\s?(\w+(?:,\s?\w+)*)', line)
    if match is None:
        raise Exception(f'{line} does not have the correct format.')
    s1, s2 = match[1], match[2]
    n2s = s2.replace(' ', '').split(',')
    return (s1, n2s)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    edge_dict: OutputT = dict()
    for line in lines:
        n1, n2s = parse_line(line)
        edge_dict[n1] = n2s

    return edge_dict


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    result = solve(kmers)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3e_1.txt'

    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    edge_dict = solve(kmers)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for seq1, seqs2 in edge_dict.items():
        out = seq1 + ' -> ' + ','.join(seqs2)
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
