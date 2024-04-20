import BioInfoToolkit.Spectrometry as Spectrometry
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Spectrometry.Spectrometry import get_prefix_masses, peptideMass


"""
https://rosalind.info/problems/spec/

The prefix spectrum of a weighted string is the collection of all its prefix weights.

    Given: A list L of n (n≤100) positive real numbers.

    Return: A protein string of length n-1 whose prefix spectrum is equal to L (if multiple solutions exist, you may output any one of them). Consult the monoisotopic mass table.
"""


def verify(result: str, solution: str, spectrum: list[float]) -> bool:
    mass_result = peptideMass(result)
    mass_solution = peptideMass(solution)

    spectrum_result = [0] + get_prefix_masses(result)
    offset = spectrum[0]
    spectrum_result = [v+offset for v in spectrum_result]

    spectrum_matches = all(abs(v1-v2) <= 0.001 for v1, v2 in zip(spectrum, spectrum_result))

    correct = len(solution) == len(result) and abs(mass_result - mass_solution) <= 0.001 and spectrum_matches

    return correct


def solve(spectrum: list[float]) -> str:
    aas = Spectrometry.getAminoacidsFromSpectrum(spectrum)
    prot = "".join(aas[0] for aas in aas)

    return prot


def load_results(path: str) -> str:
    lines = readTextFile(path)
    protein = lines[0]

    return protein


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    spectrum = [float(line) for line in lines]
    prot = solve(spectrum)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(prot, solution, spectrum)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_spec_1.txt'

    lines = readTextFile(path)
    spectrum = [float(line) for line in lines]
    prot = solve(spectrum)

    print(prot)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, prot, 'w')

    # correct = solve_and_check(path)
    # print(correct)
