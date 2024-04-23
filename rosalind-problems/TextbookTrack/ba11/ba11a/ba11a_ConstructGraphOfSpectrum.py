from BioInfoToolkit.Spectrometry.Spectrometry import getSpectrumGraph
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import re


"""
https://rosalind.info/problems/ba11a/

We represent the masses in a spectrum as a sequence Spectrum of integers  s1,…,sm in increasing order, where s1 is zero and sm is the total mass of the (unknown) peptide. We define a labeled graph Graph(Spectrum) by forming a node for each element of Spectrum, then connecting nodes si and sj by a directed edge labeled by an amino acid a if sj-si is equal to the mass of a. As we assumed when sequencing antibiotics, we do not distinguish between amino acids having the same integer masses (i.e., the pairs K/Q and I/L).

Spectrum Graph Construction
    Construct the graph of a spectrum.

        Given: A space-delimited list of integers Spectrum.

        Return: Graph(Spectrum).

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = list[tuple[int, int, str]]


def verify(result: OutputT, solution: OutputT) -> bool:
    # TODO: check if correct
    return False


def solve(spectrum: list[int]) -> OutputT:
    graph = getSpectrumGraph([float(v) for v in spectrum], tol=0.005)
    # dot = drawSpectrumGraphDot(graph)
    # dot.view()

    adj_list: list[tuple[int, int, str]] = []
    for n1, nbrs in graph.adj.items():
        w1 = graph.nodes[n1]["weight"]
        for n2, edges in nbrs.items():
            w2 = graph.nodes[n2]["weight"]
            _, eattr = next(iter(edges.items()))
            aa = eattr['label']
            adj_list.append((int(w1), int(w2), aa))

    return adj_list


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    adj_list: OutputT = []
    for line in lines:
        if not len(line) or line.isspace():
            continue

        matches = re.match(r'(\d+)->(\d+):(\w)', line)
        if matches is None:
            raise Exception('Line {line} does not have the correct format')

        n1, n2, label = int(matches[1]), int(matches[2]), matches[3]
        adj_list.append((n1, n2, label))

    return adj_list


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    spectrum = [0] + [int(val) for val in lines[0].split()]

    result = solve(spectrum)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11a_1.txt'

    lines = readTextFile(path)
    spectrum = [0] + [int(val) for val in lines[0].split()]

    adj_list = solve(spectrum)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for n1, n2, char in adj_list:
        out = f"{n1}->{n2}:{char}"
        print(out)
        writeTextFile(result_path, out, 'a')

    # correct = solve_and_check(path)
    # print(correct)
