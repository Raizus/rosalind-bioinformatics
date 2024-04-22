from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, peptideSpectrumScore
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4k/

Linear Peptide Scoring Problem
    Compute the score of a linear peptide with respect to a spectrum.

        Given: An amino acid string Peptide and a collection of integers LinearSpectrum.

        Return: The linear score of Peptide against Spectrum, LinearScore(Peptide, Spectrum).
"""


OutputT = int


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(peptide: str, spectrum: list[int]) -> OutputT:
    mass_dict = AminoacidMonoisotopicMassInteger
    score = peptideSpectrumScore(peptide, spectrum, mass_dict, False)

    return score


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    score = int(lines[0])
    return score


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    peptide = lines[0]
    spectrum = [int(m) for m in lines[1].split()]

    result = solve(peptide, spectrum)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4k_1.txt'

    lines = readTextFile(path)
    peptide = lines[0]
    spectrum = [int(m) for m in lines[1].split()]

    score = solve(peptide, spectrum)

    out = str(score)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
