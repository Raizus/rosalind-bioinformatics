from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, peptideIdentification, peptide_score
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11f/

Despite many attempts, researchers have still not devised a scoring function that reliably assigns the highest score to the biologically correct peptide, i.e., the peptide that generated the spectrum. Fortunately, although the correct peptide often does not achieve the highest score among all peptides, it typically does score highest among all peptides limited to the species's proteome. As a result, we can transition from peptide sequencing to peptide identification by limiting our search to peptides present in the proteome, which we concatenate into a single amino acid string Proteome.

Peptide Identification Problem
    Find a peptide from a proteome with maximum score against a spectrum.

        Given: A space-delimited spectral vector S and an amino acid string Proteome.

        Return: A peptide in Proteome with maximum score against S.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT, spectral_vector: list[int]) -> bool:
    score_res = peptide_score(result, spectral_vector,
                              AminoacidMonoisotopicMassInteger)
    score_sol = peptide_score(
        solution, spectral_vector, AminoacidMonoisotopicMassInteger)
    correct = score_res == score_sol
    return correct


def solve(spectral_vector: list[int], proteome: str) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger
    peptide, _ = peptideIdentification(
        spectral_vector, proteome, massDict)
    return peptide


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptide = lines[0]

    return peptide


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    spectral_vector = [0] + [int(v) for v in lines[0].split()]
    proteome = lines[1]

    result = solve(spectral_vector, proteome)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution, spectral_vector)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11f_1.txt'

    lines = readTextFile(path)
    spectral_vector = [0] + [int(v) for v in lines[0].split()]
    proteome = lines[1]

    peptide = solve(spectral_vector, proteome)

    out = peptide
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
