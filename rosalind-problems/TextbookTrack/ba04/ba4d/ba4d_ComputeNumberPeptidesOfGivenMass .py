import os
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Spectrometry.Spectrometry import countPeptides

"""
https://rosalind.info/problems/ba4d/

In “Generate the Theoretical Spectrum of a Cyclic Peptide”, we generated the theoretical spectrum of a known cyclic peptide. Although this task is relatively easy, our aim in mass spectrometry is to solve the reverse problem: we must reconstruct an unknown peptide from its experimental spectrum. We will start by assuming that a biologist is lucky enough to generate an ideal experimental spectrum Spectrum, which is one coinciding with the peptide's theoretical spectrum. Can we reconstruct a peptide whose theoretical spectrum is Spectrum?

Denote the total mass of an amino acid string Peptide as Mass(Peptide). In mass spectrometry experiments, whereas the peptide that generated a spectrum is unknown, the peptide's mass is typically known and is denoted ParentMass(Spectrum). Of course, given an ideal experimental spectrum, Mass(Peptide) is given by the largest mass in the spectrum.

A brute force approach to reconstructing a peptide from its theoretical spectrum would generate all possible peptides whose mass is equal to ParentMass(Spectrum) and then check which of these peptides has theoretical spectra matching Spectrum. However, we should be concerned about the running time of such an approach: how many peptides are there having mass equal to ParentMass(Spectrum)?

Counting Peptides with Given Mass Problem
    Compute the number of peptides of given total mass.

        Given: An integer m.

        Return: The number of linear peptides having integer mass m.
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


OutputT = int

def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(mass: int) -> OutputT:
    count = countPeptides(mass, AminoacidMonoisotopicMassInteger)
    return count


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    mass = int(lines[0])

    result = solve(mass)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4d_1.txt'

    lines = readTextFile(path)
    mass = int(lines[0])

    count = solve(mass)

    out = str(count)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
