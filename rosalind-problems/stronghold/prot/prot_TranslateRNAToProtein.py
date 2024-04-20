# https://rosalind.info/problems/dna/


from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.BioSequence import BioSequence
import os

"""
The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

    Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

    Return: The protein string encoded by s.
"""


def verify(result: str, solution: str) -> bool:
    correct = result == solution
    return correct


def solve(seq: str) -> str:
    bio_seq = BioSequence(seq, "RNA")
    prot = "".join(bio_seq.translate_seq())
    prot = "".join([aa for aa in prot if aa != '_'])
    return prot


def load_results(path: str) -> str:
    lines = readTextFile(path)
    seq = lines[0]
    return seq


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    seq = lines[0]

    result = solve(seq)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_prot_1.txt'

    lines = readTextFile(path)
    seq = lines[0]

    prot = solve(seq)
    print(prot)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, prot, 'w')
