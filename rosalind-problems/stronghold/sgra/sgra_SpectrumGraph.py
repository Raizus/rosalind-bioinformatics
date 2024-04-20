import BioInfoToolkit.Spectrometry as Spectrometry
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Spectrometry.Spectrometry import peptideMass


"""
https://rosalind.info/problems/sgra/

For a weighted alphabet A and a collection L of positive real numbers, the spectrum graph of L is a digraph constructed in the following way. First, create a node for every real number in L. Then, connect a pair of nodes with a directed edge (u,v) if v>u and v−u is equal to the weight of a single symbol in A. We may then label the edge with this symbol.

In this problem, we say that a weighted string s=s1s2⋯sn matches L if there is some increasing sequence of positive real numbers (w1,w2,…,wn+1) in L such that w(s1)=w2−w1, w(s2)=w3−w2, ..., and w(sn)=wn+1−wn.

    Given: A list L (of length at most 100) containing positive real numbers.
    
    Return: The longest protein string that matches the spectrum graph of L (if multiple solutions exist, you may output any one of them). Consult the monoisotopic mass table.
"""


def verify(result: str, solution: str) -> bool:
    correct = len(result) == len(solution) and peptideMass(result) == peptideMass(solution)

    return correct


def solve(spectrum: list[float]) -> str:
    graph = Spectrometry.getSpectrumGraph(spectrum, tol=0.001)

    prot_seqs = Spectrometry.longestSeqsFromSpectrumGraph(graph)
    prot = "".join(aas[0] for aas in prot_seqs)

    return prot


def load_results(path: str) -> str:
    lines = readTextFile(path)
    protein = lines[0]

    return protein


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    spectrum = [float(val) for val in lines]

    prot = solve(spectrum)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(prot, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_sgra_1.txt'

    lines = readTextFile(path)
    spectrum = [float(val) for val in lines]

    graph = Spectrometry.getSpectrumGraph(spectrum, tol=0.0001)

    prot_seqs = Spectrometry.longestSeqsFromSpectrumGraph(graph)
    print(prot_seqs)

    dot = Spectrometry.drawSpectrumGraphDot(graph)
    dot.view()

    prot = "".join(aas[0] for aas in prot_seqs)
    print(prot)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, prot, 'w')

    correct = solve_and_check(path)
    print(correct)
