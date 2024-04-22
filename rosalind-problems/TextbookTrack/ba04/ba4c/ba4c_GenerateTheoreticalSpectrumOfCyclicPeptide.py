from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, peptideTheoreticalSpectrum
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4c/

The workhorse of peptide sequencing is the mass spectrometer, an expensive molecular scale that shatters molecules into pieces and then weighs the resulting fragments. The mass spectrometer measures the mass of a molecule in daltons (Da); 1 Da is approximately equal to the mass of a single nuclear particle (i.e., a proton or neutron).

We will approximate the mass of a molecule by simply adding the number of protons and neutrons found in the molecule's constituent atoms, which yields the molecule's integer mass. For example, the amino acid "Gly", which has chemical formula C2H3ON, has an integer mass of 57, since 2·12 + 3·1 + 1·16 + 1·14 = 57. Yet 1 Da is not exactly equal to the mass of a proton/neutron, and we may need to account for different naturally occurring isotopes of each atom when weighing a molecule. As a result, amino acids typically have non-integer masses (e.g., "Gly" has total mass equal to approximately 57.02 Da); for simplicity, however, we will work with the integer mass table given in Figure 1.

The theoretical spectrum of a cyclic peptide Peptide, denoted Cyclospectrum(Peptide), is the collection of all of the masses of its subpeptides, in addition to the mass 0 and the mass of the entire peptide. We will assume that the theoretical spectrum can contain duplicate elements, as is the case for "NQEL" (shown in Figure 2), where "NQ" and "EL" have the same mass.

Generating Theoretical Spectrum Problem
    Generate the theoretical spectrum of a cyclic peptide.

        Given: An amino acid string Peptide.

        Return: Cyclospectrum(Peptide).
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(peptide: str) -> OutputT:
    spectrum = peptideTheoreticalSpectrum(
        peptide, AminoacidMonoisotopicMassInteger, True)
    return spectrum


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    spectrum = [int(v) for v in lines[0].split()]
    return spectrum


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    peptide = lines[0]

    result = solve(peptide)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4c_1.txt'

    lines = readTextFile(path)
    peptide = lines[0]

    spectrum = solve(peptide)

    out = ' '.join(str(v) for v in spectrum)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
