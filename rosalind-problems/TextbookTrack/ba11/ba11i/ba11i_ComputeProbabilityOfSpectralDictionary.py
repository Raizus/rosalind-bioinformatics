from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, spectralDictionaryProbability
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11i/

Note that the equation for the probability of a dictionary (introduced in “Compute the Size of a Spectral Dictionary”),

    Pr(Dictionary)=∑each peptide Peptide in Dictionary120|Peptide|,

is similar to an equation for the size of a dictionary,

    |Dictionary|=∑each peptide Peptide in Dictionary1.

This similarity suggests that we can derive a recurrence for the probability of a dictionary using arguments similar to those used to find the size of a dictionary.

So define Pr(i, t) as the sum of probabilities of all peptides with mass i for which Score(Peptide, Spectrum'i ) is equal to t. The set of peptides contributing to Pr(i, t) can be split into 20 subsets depending on their final amino acid. Each peptide Peptide ending in a specific amino acid a results in a shorter peptide Peptidea if we remove a; Peptidea has mass i - |a| and score t - si. Since the probability of Peptide is 20 times smaller than the probability of Peptidea, the contribution of Peptide to Pr(i, t) is 20 times smaller than contribution of Peptidea to Pr(i - |a|, t - si ). Therefore, Pr(i, t) can be computed as

    Pr(i,t)=∑all amino acids a120⋅Pr(i-|a|,t-si),

which differs from the recurrence for computing Size(i, t) only in the presence of the factor 1/20. We can now compute the probability of a spectral dictionary as

    Pr(Dictionarythreshold(Spectrum'))=∑t≥thresholdPr(m,t).
    

Probability of Spectral Dictionary Problem
    Find the probability of the spectral dictionary for a given spectrum and score threshold.

        Given: A spectral vector Spectrum', an integer threshold, and an integer max_score.

        Return: The probability of the dictionary Dictionarythreshold(Spectrum').

Note: Use the provided max_score for the height of your table.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = float


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = abs(result - solution) / solution <= 0.001
    return correct


def solve(spectral_vector: list[int], threshold: int, max_score: int) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger
    prob = spectralDictionaryProbability(
        spectral_vector, threshold, max_score, massDict)
    return prob


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    prob = float(lines[0])

    return prob


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)

    spectral_vector = [0] + [int(v) for v in lines[0].split()]
    threshold = int(lines[1])
    max_score = int(lines[2])

    result = solve(spectral_vector, threshold, max_score)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11i_1.txt'

    lines = readTextFile(path)

    spectral_vector = [0] + [int(v) for v in lines[0].split()]
    threshold = int(lines[1])
    max_score = int(lines[2])

    prob = solve(spectral_vector, threshold, max_score)

    out = str(prob)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
