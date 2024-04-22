from collections import Counter
from BioInfoToolkit.Sequences.GenomeAssembly import OverlapGraph
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import re

"""
https://rosalind.info/problems/ba3c/

String Spelled by a Genome Path Problem
    Find the string spelled by a genome path.

        Given: A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1 symbols of Patterni are equal to the first k - 1 symbols of Patterni+1 for i from 1 to n-1.

        Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.
"""

OutputT = list[tuple[str, str]]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(kmers: list[str]) -> OutputT:
    k = len(kmers[0])

    overlapGraph = OverlapGraph(kmers, k-1)
    # dot = overlapGraph.draw_dot()
    # dot.view()
    edge_list: OutputT = []
    for edge in overlapGraph.graph.edges:
        n1, n2 = edge
        seq1 = overlapGraph.graph.nodes[n1]['seq']
        seq2 = overlapGraph.graph.nodes[n2]['seq']
        edge_list.append((seq1, seq2))

    return edge_list


def parse_line(line: str) -> tuple[str,str]:
    match = re.match(r'(\w+)\s?->\s?(\w+)', line)
    if match is None:
        raise Exception(f'{line} does not have the correct format.')
    s1, s2 = match[1], match[2]
    return (s1,s2)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    edge_list: OutputT = []
    for line in lines:
        edge = parse_line(line)
        edge_list.append(edge)

    return edge_list


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
    path = f'{cwd}/rosalind_ba3c_1.txt'

    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    edge_list = solve(kmers)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for seq1, seq2 in edge_list:
        out = seq1 + ' -> ' + seq2
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)

