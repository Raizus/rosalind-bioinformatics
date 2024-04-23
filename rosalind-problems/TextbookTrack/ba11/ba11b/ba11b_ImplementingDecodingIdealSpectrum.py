from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, SpectrumGraph, peptideTheoreticalSpectrum
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11b/

Given an amino acid string Peptide, its ideal spectrum, denoted IdealSpectrum(Peptide), is the collection of integer masses of all its prefixes and suffixes. Note that an ideal spectrum may have repeated masses; for example, IdealSpectrum(GPG) = {0, 57, 57, 154, 154, 211}. We say that an amino acid string Peptide explains a collection of integers Spectrum if IdealSpectrum(Peptide) = Spectrum.

The following pseudocode finds a peptide explaining a given spectrum.
    DecodingIdealSpectrum(Spectrum)
        construct Graph(Spectrum)
        for each path Path from source to sink in Graph(Spectrum)
            Peptide ← the amino acid string spelled by the edge labels of Path
            if IdealSpectrum(Peptide) = Spectrum
                    return Peptide


Decoding an Ideal Spectrum Problem
    Reconstruct a peptide from its ideal spectrum.

        Given: A space-delimited list of integers, Spectrum.

        Return: An amino acid string with an ideal spectrum that matches Spectrum.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    mass_dict = AminoacidMonoisotopicMassInteger
    spectrum_res = peptideTheoreticalSpectrum(result, mass_dict, False)
    spectrum_sol = peptideTheoreticalSpectrum(solution, mass_dict, False)
    correct = spectrum_res == spectrum_sol
    return correct


def solve(spectrum: list[int]) -> OutputT:
    spectrumGrap = SpectrumGraph([float(v) for v in spectrum], tol=0.005)
    # dot = spectrumGrap.draw_dot()
    # dot.view()

    peptide = spectrumGrap.decode_ideal_spectrum()

    return peptide


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptide = lines[0]

    return peptide


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
    path = f'{cwd}/rosalind_ba11b_1.txt'

    lines = readTextFile(path)
    spectrum = [0] + [int(val) for val in lines[0].split()]

    peptide = solve(spectrum)

    out = peptide
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
