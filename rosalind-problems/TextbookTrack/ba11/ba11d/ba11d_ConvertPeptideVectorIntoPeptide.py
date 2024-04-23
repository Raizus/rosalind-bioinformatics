from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, get_peptide_from_vector, get_peptide_vector
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11d/

We can also invert the conversion of a peptide into a peptide vector, so that the terms "peptide" and "peptide vector" can be used interchangeably.

Converting Peptide Vector into Peptide Problem
    Convert a binary vector into a peptide.

        Given: A space-delimited binary vector P.

        Return: A peptide whose binary peptide vector matches P. For masses with more than one amino acid, any choice may be used.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    peptide_vector_res = get_peptide_vector(
        result, AminoacidMonoisotopicMassInteger)
    peptide_vector_sol = get_peptide_vector(
        solution, AminoacidMonoisotopicMassInteger)
    correct = peptide_vector_res == peptide_vector_sol
    return correct


def solve(peptide_vector: list[int]) -> OutputT:
    peptide = get_peptide_from_vector(peptide_vector)
    return peptide


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptide = lines[0]

    return peptide


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    peptide_vetor = [int(v) for v in lines[0].split()]

    result = solve(peptide_vetor)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11d_1.txt'

    lines = readTextFile(path)
    peptide_vetor = [int(v) for v in lines[0].split()]

    peptide = solve(peptide_vetor)

    out = peptide
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
