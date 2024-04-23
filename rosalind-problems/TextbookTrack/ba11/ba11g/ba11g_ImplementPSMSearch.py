from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, PSMSearch, peptide_score
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11g/

Like peptide sequencing algorithms, peptide identification algorithms may return an erroneous peptide, particularly if the score of the highest-scoring peptide found in the proteome is much lower than the score of the highest-scoring peptide over all peptides. For this reason, biologists usually establish a score threshold and only pay attention to a solution of the Peptide Identification Problem if its score is at least equal to the threshold.

Given a set of spectral vectors SpectralVectors, an amino acid string Proteome, and a score threshold threshold, we will solve the Peptide Identification Problem for each vector Spectrum' in SpectralVectors and identify a peptide Peptide having maximum score for this spectral vector over all peptides in Proteome (ties are broken arbitrarily). If Score(Peptide, Spectrum) is greater than or equal to threshold, then we conclude that Peptide is present in the sample and call the pair (Peptide, Spectrum') a Peptide- Spectrum Match (PSM). The resulting collection of PSMs for SpectralVectors is denoted PSMthreshold(Proteome, SpectralVectors).

The following pseudocode identifies all PSMs scoring above a threshold for a set of spectra and a proteome, using an algorithm that you just implemented to solve the Peptide identification Problem, which we call PeptideIdentification.

    PSMSearch(SpectralVectors, Proteome, threshold).
        PSMSet ← an empty set
        for each vector Spectrum' in SpectralVectors
            Peptide ← PeptideIdentification(Spectrum', Proteome)
            if Score(Peptide, Spectrum) ≥ threshold
                add the PSM (Peptide, Spectrum') to PSMSet
        return PSMSet

        
PSM Search Problem
    Identify Peptide-Spectrum Matches by matching spectra against a proteome.

        Given: A set of space-delimited spectral vectors SpectralVectors, an amino acid string Proteome, and a score threshold T.

        Return: All unique Peptide-Spectrum Matches scoring at least as high as T.

Note: For this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X: 4, Z: 5.
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT, spectral_vectors: list[list[int]], 
           proteome: str, threshold: int) -> bool:

    result_above_threshold = all(any(peptide_score(peptide, spectrum, AminoacidMonoisotopicMassInteger) >= threshold for spectrum in spectral_vectors) for peptide in result)

    peps_in_proteome = all(proteome.find(peptide) != -1 for peptide in result )

    correct = len(solution) == len(
        result) and peps_in_proteome and result_above_threshold and set(result) == set(solution)
    return correct


def solve(spectral_vectors: list[list[int]], proteome: str, threshold: int) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger
    PSMSet = PSMSearch(spectral_vectors, proteome, threshold, massDict)
    return list(PSMSet)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptides = [line for line in lines if len(line) and not line.isspace()]

    return peptides


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)

    spectral_vectors: list[list[int]] = []
    last_i = 0
    for i, line in enumerate(lines):
        res = line.split()
        if len(res) > 1:
            spectral_vector = [0] + [int(v) for v in res]
            spectral_vectors.append(spectral_vector)
        else:
            last_i = i
            break
    proteome = lines[last_i]
    threshold = int(lines[last_i+1])

    result = solve(spectral_vectors, proteome, threshold)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution, spectral_vectors, proteome, threshold)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11g_1.txt'

    lines = readTextFile(path)

    spectral_vectors: list[list[int]] = []
    last_i = 0
    for i, line in enumerate(lines):
        res = line.split()
        if len(res) > 1:
            spectral_vector = [0] + [int(v) for v in res]
            spectral_vectors.append(spectral_vector)
        else:
            last_i = i
            break
    proteome = lines[last_i]
    threshold = int(lines[last_i+1])

    peptides = solve(spectral_vectors, proteome, threshold)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for peptide in peptides:
        print(peptide)
        writeTextFile(result_path, peptide, 'a')

    correct = solve_and_check(path)
    print(correct)
