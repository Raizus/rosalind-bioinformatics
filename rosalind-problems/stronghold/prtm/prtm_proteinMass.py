from BioInfoToolkit.Spectrometry.Spectrometry import peptideMass
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
In a weighted alphabet, every symbol is assigned a positive real number called a weight. A string formed from a weighted alphabet is called a weighted string, and its weight is equal to the sum of the weights of its symbols.

The standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic mass of the corresponding amino acid.

    Given: A protein string P of length at most 1000 aa.

    Return: The total weight of P. Consult the monoisotopic mass table.
"""


def verify(result: float, solution: float) -> bool:
    error = abs(result - solution)
    correct = error <= 0.001
    return correct


def solve(peptide: str) -> float:
    mass = peptideMass(peptide)
    return mass


def load_results(path: str) -> float:
    lines = readTextFile(path)
    mass = float(lines[0])
    return mass


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    peptide = lines[0]

    mass = solve(peptide)
    solution = load_results(solution_path)

    correct = verify(mass, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_prtm_1.txt'

    lines = readTextFile(path)
    peptide = lines[0]

    mass = solve(peptide)

    print(round(mass, 4))

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(round(mass, 4)), 'w')

