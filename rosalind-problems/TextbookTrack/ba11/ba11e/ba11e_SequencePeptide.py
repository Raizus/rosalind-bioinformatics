from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, SpectralGraph, peptide_score
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11e/

Given a spectral vector Spectrum', our goal is to find a peptide Peptide maximizing Score(Peptide, Spectrum'). Since the mass of a peptide and the parent mass of the spectrum that it generates should be the same, a peptide vector should have the same length as the spectral vector under consideration. We will therefore define the score between a peptide vector and a spectral vector of different length as −∞.

Given a spectral vector Spectrum' = (s1, . . . , sm), our goal is to find a peptide with maximum score against this spectral vector. To do so, we will construct a DAG on m + 1 nodes, labeled with the integers from 0 (source) to m (sink), and then connect node i to node j by a directed edge if j − i is equal to the mass of an amino acid. We will further assign weight si to node i (for 1 ≤ i ≤ m) and assign weight zero to node 0.

Any path connecting source to sink in this DAG corresponds to an amino acid string Peptide, and the total weight of nodes on this path is equal to Score(Peptide', Spectrum'). We have therefore reduced the Peptide Sequencing Problem to the problem of finding a maximum-weight path from source to sink in a node-weighted DAG.

Peptide Sequencing Problem
    Given a spectral vector S, find a peptide vector with maximum score against S.

        Given: A space-delimited spectral vector S.

        Return: A peptide with maximum score against S. For masses with more than one amino acid, any choice may be used.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT, spectral_vector: list[int]) -> bool:
    score_res = peptide_score(result, spectral_vector, AminoacidMonoisotopicMassInteger)
    score_sol = peptide_score(solution, spectral_vector, AminoacidMonoisotopicMassInteger)
    correct = score_res == score_sol
    return correct


def solve(spectral_vector: list[int]) -> OutputT:
    spectalGraph = SpectralGraph(
        spectral_vector, 0.005, AminoacidMonoisotopicMassInteger)
    # dot = spectalGraph.draw_dot()
    # dot.view('spectral_graph.svg')

    peptide = spectalGraph.peptide_sequence()
    return peptide


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptide = lines[0]

    return peptide


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    spectral_vector = [0] + [int(v) for v in lines[0].split()]

    result = solve(spectral_vector)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution, spectral_vector)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11e_0.txt'

    lines = readTextFile(path)
    spectral_vector = [0] + [int(v) for v in lines[0].split()]

    peptide = solve(spectral_vector)

    out = peptide
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
