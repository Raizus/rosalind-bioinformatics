# https://rosalind.info/problems/full/

import math
from typing import Iterable
import BioInfoToolkit.Spectrometry as Spectrometry
from BioInfoToolkit.IO import readTextFile
import networkx as nx

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Spectrometry.Spectrometry import peptideMass


"""
https://rosalind.info/problems/full/

Say that we have a string s containing t as an internal substring, so that there exist nonempty substrings s1 and s2 of s such that s can be written as s1ts2. A t-prefix contains all of s1 and none of s2; likewise, a t-suffix contains all of s2 and none of s1.

    Given: A list L containing 2n+3 positive real numbers (n≤100). The first number in L is the parent mass of a peptide P, and all other numbers represent the masses of some b-ions and y-ions of P (in no particular order). You may assume that if the mass of a b-ion is present, then so is that of its complementary y-ion, and vice-versa.

    Return: A protein string t of length n for which there exist two positive real numbers w1 and w2 such that for every prefix p and suffix s of t, each of w(p)+w1 and w(s)+w2 is equal to an element of L. (In other words, there exists a protein string whose t-prefix and t-suffix weights correspond to the non-parent mass values of L.) If multiple solutions exist, you may output any one.
"""


def findPathsOfLenghtN(g: nx.MultiDiGraph, node: int, n: int) -> Iterable[list[tuple[int, int, int]]]:
    if n == 0:
        yield []
    for node2, edges in g.adj[node].items():
        for edgeIdx, _ in edges.items():
            edge: tuple[int, int, int] = (node, node2, edgeIdx)
            for path in findPathsOfLenghtN(g, node2, n-1):
                yield [edge] + path


def findPathLowestError(g: nx.MultiDiGraph, node: int, n: int):
    error = math.inf
    shortestPath = None
    for path in findPathsOfLenghtN(graph, node, n):
        _error = sum(g.adj[u][v][e]['error'] for (u, v, e) in path)
        if _error < error:
            error = _error
            shortestPath = path
    return shortestPath


def pathToSeq(g: nx.MultiDiGraph, path: list[tuple[int, int, int]]) -> str:
    seq = "".join(g.adj[u][v][e]['label'] for (u, v, e) in path)
    return seq



def verify(result: str, solution: str) -> bool:
    mass_result = peptideMass(result)
    mass_solution = peptideMass(solution)

    correct = len(solution) == len(result) and abs(
        mass_result - mass_solution) <= 0.001

    return correct


def solve(peptide_mass: float, spectrum: list[float]) -> str:
    n = len(spectrum)//2 - 1
    graph = Spectrometry.getSpectrumGraph(spectrum, tol=0.0001)
    Spectrometry.longestSeqsFromSpectrumGraph(graph)
    
    graph_path = findPathLowestError(graph, 0, n)
    if graph_path:
        prot = pathToSeq(graph, graph_path)
        return prot

    raise Exception("No path found in spectrum graph")

def load_results(path: str) -> str:
    lines = readTextFile(path)
    protein = lines[0]

    return protein


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    masses = [float(line) for line in lines]
    peptide_mass = masses[0]
    full_spectrum = sorted(masses[1:])

    prot = solve(peptide_mass, full_spectrum)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(prot, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_full_1.txt'

    lines = readTextFile(path)
    masses = [float(line) for line in lines]
    peptide_mass = masses[0]
    full_spectrum = sorted(masses[1:])

    n = len(full_spectrum)//2 - 1
    graph = Spectrometry.getSpectrumGraph(full_spectrum, tol=0.0001)
    Spectrometry.longestSeqsFromSpectrumGraph(graph)
    
    # dot = Spectrometry.drawSpectrumGraphDot(graph)
    # dot.view()

    graph_path = findPathLowestError(graph, 0, n)
    if graph_path:
        seq = pathToSeq(graph, graph_path)
        print(seq)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, seq, 'w')

    correct = solve_and_check(path)
    print(correct)