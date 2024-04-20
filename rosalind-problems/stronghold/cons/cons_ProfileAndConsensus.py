from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

from BioInfoToolkit.Sequences.Profiling import SequencesProfile
from BioInfoToolkit.Sequences.structures import NUCLEOTIDE_BASE
import re

"""
A matrix is a rectangular table of values divided into rows and columns. An m x n matrix has m rows and n columns. Given a matrix A, we write Ai,j to indicate the value found at the intersection of row i and column j.

Say that we have a collection of DNA strings, all having the same length n. Their profile matrix is a 4 x n matrix P in which P1,j represents the number of times that 'A' occurs in the jth position of one of the strings, P2,j represents the number of times that C occurs in the jth position, and so on (see below).

A consensus string c is a string of length n formed from our collection by taking the most common symbol at each position; the jth symbol of c therefore corresponds to the symbol having the maximum value in the j-th column of the profile matrix. Of course, there may be more than one most common symbol, leading to multiple possible consensus strings.

	        A T C C A G C T
	        G G G C A A C T
	        A T G G A T C T
DNA Strings	A A G C A A C C
	        T T G G A A C T
            A T G C C A T T
            A T G G C A C T

            A   5 1 0 0 5 5 0 0
Profile	    C   0 0 1 4 2 0 6 1
        	G   1 1 6 3 0 1 0 0
        	T   1 5 0 0 0 1 1 6
Consensus	    A T G C A A C T

    Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

    Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)
"""


def verify(result: tuple[str, SequencesProfile], solution: tuple[str, SequencesProfile]) -> bool:
    result_profile = result[1]
    result_consensus = result[0]
    solution_profile = solution[1]

    correct_profile = all(rval == sval 
                          for res_line, sol_line in zip(result_profile.count_matrix, solution_profile.count_matrix) 
                          for rval, sval in zip(res_line, sol_line))

    if not correct_profile:
        return False
    
    is_consensus = result_profile.is_consensus(result_consensus)
    return is_consensus


def solve(fasta_dict: dict[str, str]) -> tuple[SequencesProfile, str]:
    seqs = list(fasta_dict.values())
    profile = SequencesProfile(seqs, NUCLEOTIDE_BASE['DNA'])
    consensus = profile.consensus()
    return profile, consensus


def load_results(path: str) -> tuple[str, SequencesProfile]:
    lines = readTextFile(path)
    consensus = lines[0]
    profile_lines = [line for line in lines[1:] if len(line)]

    alphabet: list[str] = []
    profile_mat: list[list[int]] = []
    for line in profile_lines:
        matches = re.match(r'(\w):((?:\s\d+)+)', line)
        if not matches:
            raise Exception('Result file does not have the correct format.')
        alphabet.append(matches[1])
        counts = [int(v) for v in matches[2].split()]
        profile_mat.append(counts)
    profile = SequencesProfile.from_profile_matrix(alphabet, profile_mat)

    return consensus, profile


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    profile, consensus = solve(fasta_dict)

    solution = load_results(solution_path)

    correct = verify((consensus, profile), solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_cons_1.txt'

    fasta_dict = read_FASTA(path)
    profile, consensus = solve(fasta_dict)

    print(consensus)
    print(profile)

    result_path = result_path_from_input_path(path)
    out = consensus + '\n' + str(profile)
    writeTextFile(result_path, out, 'w')

