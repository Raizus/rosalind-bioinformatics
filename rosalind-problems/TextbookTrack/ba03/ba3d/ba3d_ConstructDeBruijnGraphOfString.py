from collections import Counter
from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromString
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import re

"""
https://rosalind.info/problems/ba3d/

Given a genome Text, PathGraphk(Text) is the path consisting of |Text| - k + 1 edges, where the i-th edge of this path is labeled by the i-th k-mer in Text and the i-th node of the path is labeled by the i-th (k - 1)-mer in Text. The de Bruijn graph DeBruijnk(Text) is formed by gluing identically labeled nodes in PathGraphk(Text).

De Bruijn Graph from a String Problem
    Construct the de Bruijn graph of a string.

        Given: An integer k and a string Text.

        Return: DeBruijnk(Text), in the form of an adjacency list.
"""

OutputT = dict[str, list[str]]


def verify(result: OutputT, solution: OutputT) -> bool:
    if set(result.keys()) != set(solution.keys()):
        return False
    
    for key1, val1 in result.items():
        if Counter(solution[key1]) != Counter(val1):
            return False
    return True


def solve(sequence: str, k: int) -> OutputT:
    _g = deBruijnMultiGraphFromString(sequence, k)
    g = DeBruijnMultiGraph(_g, k)
    # dot = g.draw_dot()
    # dot.view()

    # string = g.reconstructStringFromEulerianPath()
    # print(string)

    sorted_nodes = sorted([nodeId for nodeId in g.graph.nodes])
    edge_dict: OutputT = dict()
    for node1 in sorted_nodes:
        adjacency_dict = g.graph.adj[node1]
        adj_seqs: list[str] = list(adjacency_dict.keys())

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
    k = int(lines[0])
    text = lines[1]

    result = solve(text, k)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3d_1.txt'

    lines = readTextFile(path)
    k = int(lines[0])
    text = lines[1]

    edge_dict = solve(text, k)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for seq1, seqs2 in edge_dict.items():
        out = seq1 + ' -> ' + ','.join(seqs2)
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
