from collections import Counter
from BioInfoToolkit.Sequences.StringUtils import kmer_gen
from BioInfoToolkit.Spectrometry.Spectrometry import cyclopeptideSequencing
import os
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba4e/

In “Compute the Number of Peptides of Given Total Mass”, we first encountered the problem of reconstructing a cyclic peptide from its theoretical spectrum; this problem is called the Cyclopeptide Sequencing Problem and is given below. It is solved by the following algorithm.

    CYCLOPEPTIDESEQUENCING(Spectrum)
        Peptides ← a set containing only the empty peptide
        while Peptides is nonempty
            Peptides ← Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides

Cyclopeptide Sequencing Problem
    Given an ideal experimental spectrum, find a cyclic peptide whose theoretical spectrum matches the experimental spectrum.

        Given: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.

        Return: Every amino acid string Peptide such that Cyclospectrum(Peptide) = Spectrum (if such a string exists).
"""

# have to remove duplicate masses for this problem, which is stupid
AminoacidMonoisotopicMassInteger = {
    "A":   71,
    "C":   103,
    "D":   115,
    "E":   129,
    "F":   147,
    "G":   57,
    "H":   137,
    "I":   113,
    "K":   128,
    "M":   131,
    "N":   114,
    "P":   97,
    "R":   156,
    "S":   87,
    "T":   101,
    "V":   99,
    "W":   186,
    "Y":   163}



OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(spectrum: list[int]) -> OutputT:
    mass_dict = AminoacidMonoisotopicMassInteger
    peptide = cyclopeptideSequencing(spectrum, mass_dict)
    
    res: list[str] = []
    if peptide:
        for kmer in kmer_gen(peptide, len(peptide), True):
            text = '-'.join(str(mass_dict[aa]) for aa in kmer)
            res.append(text)
        for kmer in kmer_gen(peptide[::-1], len(peptide), True):
            text = '-'.join(str(mass_dict[aa]) for aa in kmer)
            res.append(text)

    return res


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptides = lines[0].split()
    return peptides


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    spectrum = [int(m) for m in lines[0].split()]

    result = solve(spectrum)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4e_1.txt'

    lines = readTextFile(path)
    spectrum = [int(m) for m in lines[0].split()]

    peptides = solve(spectrum)

    out = ' '.join(peptides)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
