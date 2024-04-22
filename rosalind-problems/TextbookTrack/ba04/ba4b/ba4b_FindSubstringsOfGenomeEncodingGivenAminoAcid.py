from collections import Counter
from BioInfoToolkit.Sequences.BioSequence import BioSequence
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4b/

There are three different ways to divide a DNA string into codons for translation, one starting at each of the first three starting positions of the string. These different ways of dividing a DNA string into codons are called reading frames. Since DNA is double-stranded, a genome has six reading frames (three on each strand), as shown in Figure 1.

We say that a DNA string Pattern encodes an amino acid string Peptide if the RNA string transcribed from either Pattern or its reverse complement Pattern translates into Peptide.

Peptide Encoding Problem
    Find substrings of a genome encoding a given amino acid sequence.

        Given: A DNA string Text and an amino acid string Peptide.

        Return: All substrings of Text encoding Peptide (if any such substrings exist).
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(sequence: str, peptide: str) -> OutputT:
    bioSeq = BioSequence(sequence, "DNA")
    matches = bioSeq.peptideDecoding(peptide)
    return matches


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    sequences = [line for line in lines if len(line) and not line.isspace()]
    return sequences


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    dna_string = lines[0]
    peptide = lines[1]

    result = solve(dna_string, peptide)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4b_1.txt'

    lines = readTextFile(path)
    dna_string = lines[0]
    peptide = lines[1]

    encoding_strings = solve(dna_string, peptide)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for string in encoding_strings:
        print(string)
        writeTextFile(result_path, string, 'a')

    correct = solve_and_check(path)
    print(correct)

