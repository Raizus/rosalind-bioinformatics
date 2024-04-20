# https://rosalind.info/problems/dna/


from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.BioSequence import BioSequence
from BioInfoToolkit.Sequences.structures import DNA_NUCLEOTIDES
import os

"""
A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

    Given: A DNA string s of length at most 1000 nt.

    Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.
"""

def verify(result: list[int], solution: list[int]) -> bool:
    correct = all(r == s for r, s in zip(result, solution))
    return correct


def solve(seq: str) -> list[int]:
    bioSeq = BioSequence(seq)
    freq = bioSeq.nucleotide_frequency()
    vals = [freq[nuc] for nuc in DNA_NUCLEOTIDES]

    return vals


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    result = [int(v) for v in lines[0].split()]
    return result


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    result = solve(seq)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct

if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_dna_1.txt'

    lines = readTextFile(path)
    seq = lines[0]

    result = solve(seq)
    out = " ".join(str(val) for val in result)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')
    # solve_and_check(path)
