from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, spectralDictionarySize
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/ba11h/

Imagine that we have an algorithm generating the set of all peptides scoring at least threshold against a spectral vector Spectrum'. We will henceforth call this set of high-scoring peptides a spectral dictionary, denoted Dictionarythreshold(Spectrum'). For a PSM (Peptide, Spectrum'), we will use the term PSM dictionary, denoted Dictionary(Peptide, Spectrum'), to refer to the spectral dictionary DictionaryScore(Peptide, Spectrum')(Spectrum).

Say we have a spectral dictionary consisting of a single amino acid string Peptide, which we will attempt to match against a randomly generated string DecoyProteome of length n. Because DecoyProteome was randomly generated, the probability that Peptide matches the string beginning at a given position of DecoyProteome is 1/20|Peptide|. We call this expression the probability of Peptide.

If the spectral dictionary consists of more than one peptide, then we define the probability of the dictionary as

    ∑each peptide Peptide in Dictionary120|Peptide|.

We will first compute the number of peptides in a spectral dictionary, since this simpler problem will provide insights on how to compute the probability of a spectral dictionary.

We will use dynamic programming to find the size of a spectral Dictionary problem. Given a spectral vector Spectrum' = (s1, . . . , sm), we define its i-prefix (for i between 1 and m) as Spectrum'i = (s1, . . . , si ) and introduce a variable Size(i, t) as the number of peptides Peptide of mass i such that Score(Peptide, Spectrum'i) is equal to t.

The key to establishing a recurrence relation for computing Size(i, t) is to realize that the set of peptides contributing to Size(i, t) can be split into 20 subsets depending on their final amino acid a. Each peptide ending in a specific amino acid a results in a shorter peptide with mass i - |a| and score t - si   if we remove a from the peptide (here, |a| denotes the mass of a). Thus,

    Size(i,t)=∑all amino acids aSize(i-|a|,t-si).

Since there is a single “empty” peptide of length zero, we initialize Size(0, 0) = 1. We also define Size(0, t) = 0 for all possible scores t, and set Size(i, t) = 0 for negative values of i. Using the above recurrence, we can compute the size of a spectral dictionary of Spectrum' = (s1, . . . , sm) as

    |Dictionarythreshold(Spectrum)|=∑t≥thresholdSize(m,t).


Size of Spectral Dictionary Problem
    Find the size of the spectral dictionary for a given spectrum and score threshold.

        Given: A spectral vector Spectrum', an integer threshold, and an integer max_score.

        Return: The size of the dictionary Dictionarythreshold(Spectrum').

Note: Use the provided max_score for the height of your table. Your answer should be the number of peptides whose score is at least threshold and at most max_score.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = int


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(spectral_vector: list[int], threshold: int, max_score: int) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger
    size = spectralDictionarySize(
        spectral_vector, threshold, max_score, massDict)
    return size


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    size = int(lines[0])

    return size


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
    path = f'{cwd}/rosalind_ba11h_1.txt'

    lines = readTextFile(path)

    spectral_vector = [0] + [int(v) for v in lines[0].split()]
    threshold = int(lines[1])
    max_score = int(lines[2])

    size = solve(spectral_vector, threshold, max_score)

    out = str(size)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
