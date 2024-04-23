from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, get_peptide_vector
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11c/

Given an amino acid string Peptide = a1 . . . an of length n, we will represent its prefix masses using a binary peptide vector Peptide' with mass(Peptide) coordinates. This vector contains a 1 at each of the n prefix coordinates

    mass(a1), mass(a1 a2), . . . , mass(a1 a2 . . . an ) ,

and it contains a 0 in each of the remaining noise coordinates. The toy peptide XZZXX, whose prefix masses are 4, 9, 14, 18, and 22, corresponds to the peptide vector (0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1) of length 22.

Converting Peptide into Peptide Vector Problem
    Convert a peptide into a binary peptide vector.

        Given: A peptide P.

        Return: The peptide vector of P.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(peptide: str) -> OutputT:
    peptide_vector = get_peptide_vector(
        peptide, AminoacidMonoisotopicMassInteger)

    return peptide_vector


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptide_vector = [int(v) for v in lines[0].split()]

    return peptide_vector


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    peptide = lines[0]

    result = solve(peptide)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11c_1.txt'

    lines = readTextFile(path)
    peptide = lines[0]

    peptide_vector = solve(peptide)

    out = ' '.join(str(v) for v in peptide_vector)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
