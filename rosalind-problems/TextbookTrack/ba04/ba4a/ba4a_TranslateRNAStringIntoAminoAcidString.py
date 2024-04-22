from BioInfoToolkit.Sequences.BioSequence import BioSequence
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4a/

Much like replication, the chemical machinery underlying transcription and translation is fascinating, but from a computational perspective, both processes are straightforward. Transcription simply transforms a DNA string into an RNA string by replacing all occurrences of "T" with "U". The resulting strand of RNA is translated into an amino acid sequence via the genetic code; this process converts each 3-mer of RNA, called a codon, into one of 20 amino acids.

As illustrated in Figure 1, each of the 64 RNA codons encodes its own amino acid (some codons encode the same amino acid), with the exception of three stop codons that do not translate into amino acids and serve to halt translation. For example, the DNA string "TATACGAAA" transcribes into the RNA string "UAUACGAAA", which in turn translates into the amino acid string "Tyr-Thr-Lys".

The following problem asks you to find the translation of an RNA string into an amino acid string.

Protein Translation Problem
    Translate an RNA string into an amino acid string.

        Given: An RNA string Pattern.

        Return: The translation of Pattern into an amino acid string Peptide.
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(sequence: str) -> OutputT:
    bioSeq = BioSequence(sequence, "RNA")
    peptide = ''.join(bioSeq.translate_seq())
    peptide = peptide[:-1]
    return peptide


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    protein = lines[0]
    return protein


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    sequence = lines[0]

    result = solve(sequence)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4a_2.txt'

    lines = readTextFile(path)
    sequence = lines[0]

    protein = solve(sequence)

    out = protein
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

