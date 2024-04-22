from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, peptideTheoreticalSpectrum
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4j/

Given an amino acid string Peptide, we will begin by assuming that it represents a linear peptide. Our approach to generating its theoretical spectrum is based on the assumption that the mass of any subpeptide is equal to the difference between the masses of two prefixes of Peptide. We can compute an array PrefixMass storing the masses of each prefix of Peptide in increasing order, e.g., for Peptide = NQEL, PrefixMass = (0, 114, 242, 371, 484). Then, the mass of the subpeptide of Peptide beginning at position i + 1 and ending at position j can be computed as PrefixMass(j) - PrefixMass(i). For example, when Peptide = NQEL,

Mass(QE) = PrefixMass(3) - PrefixMass(1) = 371 - 114 = 257.

The pseudocode shown on the next step implements this idea. It also represents the alphabet of 20 amino acids and their integer masses as a pair of 20-element arrays AminoAcid and AminoAcidMass, corresponding to the top and bottom rows of the following integer mass table, respectively.

    LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)
        PrefixMass(0) ← 0
        for i ← 1 to |Peptide|
            for j ← 1 to 20
                if AminoAcid(j) =  i-th amino acid in Peptide
                    PrefixMass(i) ← PrefixMass(i - 1) + AminoAcidMass(j)
        LinearSpectrum ← a list consisting of the single integer 0
        for i ← 0 to |Peptide| - 1
            for j ← i + 1 to |Peptide|
                add PrefixMass(j) - PrefixMass(i) to LinearSpectrum
        return sorted list LinearSpectrum

Linear Spectrum Problem
    Generate the ideal linear spectrum of a peptide.

        Given: An amino acid string Peptide.

        Return: The linear spectrum of Peptide.
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(peptide: str) -> OutputT:
    mass_dict = AminoacidMonoisotopicMassInteger
    spectrum = peptideTheoreticalSpectrum(peptide, mass_dict, False)
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
    path = f'{cwd}/rosalind_ba4j_1.txt'

    lines = readTextFile(path)
    peptide = lines[0]

    spectrum = solve(peptide)

    out = ' '.join(str(v) for v in spectrum)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
