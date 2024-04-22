from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, peptideSpectrumScore
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4f/

To generalize the Cyclopeptide Sequencing Problem from “Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum” to handle noisy spectra, we need to relax the requirement that a candidate peptide's theoretical spectrum must match the experimental spectrum exactly, and instead incorporate a scoring function that will select the peptide whose theoretical spectrum matches the given experimental spectrum the most closely. Given a cyclic peptide Peptide and a spectrum Spectrum, we define Score(Peptide, Spectrum) as the number of masses shared between Cyclospectrum(Peptide) and Spectrum. Recalling our example above, if

    >Spectrum = {0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484},

then Score(NQEL, Spectrum) = 11.

The scoring function should take into account the multiplicities of shared masses, i.e., how many times they occur in each spectrum. For example, suppose that Spectrum is the theoretical spectrum of NQEL; for this spectrum, mass 242 has multiplicity 2. If 242 has multiplicity 1 in the theoretical spectrum of Peptide, then 242 contributes 1 to Score(Peptide, Spectrum). If 242 has larger multiplicity in the theoretical spectrum of Peptide, then 242 contributes 2 to Score(Peptide, Spectrum).

Cyclic Peptide Scoring Problem
    Compute the score of a cyclic peptide against a spectrum.

        Given: An amino acid string Peptide and a collection of integers Spectrum.

        Return: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
"""


OutputT = int


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(peptide: str, spectrum: list[int]) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger
    score = peptideSpectrumScore(peptide, spectrum, massDict, True)

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
    path = f'{cwd}/rosalind_ba4f_1.txt'

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
