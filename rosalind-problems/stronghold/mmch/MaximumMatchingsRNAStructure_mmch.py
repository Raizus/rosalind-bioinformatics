
from collections import Counter
import math
import os
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
# from BioInfoToolkit.Sequences.SequenceUtils import bondingGraph, drawBondingGraph
# from BioInfoToolkit.Sequences.structures import RNA_COMPLEMENT

from BioInfoToolkit.IO import read_FASTA

"""
https://rosalind.info/problems/mmch/

The graph theoretical analogue of the quandary stated in the introduction above is that if we have an RNA string s that does not have the same number of occurrences of 'C' as 'G' and the same number of occurrences of 'A' as 'U', then the bonding graph of s cannot possibly possess a perfect matching among its basepair edges. For example, see Figure 1; in fact, most bonding graphs will not contain a perfect matching.

In light of this fact, we define a maximum matching in a graph as a matching containing as many edges as possible. See Figure 2 for three maximum matchings in graphs.

A maximum matching of basepair edges will correspond to a way of forming as many base pairs as possible in an RNA string, as shown in Figure 3.

    Given: An RNA string s of length at most 100.

    Return: The total possible number of maximum matchings of basepair edges in the bonding graph of s.
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq: str) -> int:
    counts = Counter(seq)
    minAU = min(counts['A'], counts['U'])
    maxAU = max(counts['A'], counts['U'])
    minGC = min(counts['G'], counts['C'])
    maxGC = max(counts['G'], counts['C'])

    # The number of perfect matches will be equal to the all the possible ways to match all the A's to all the U's times the number of possible ways to match all the G's to all the C's 
    numMaximumMatching = math.perm(maxAU, minAU) * math.perm(maxGC, minGC)
    return numMaximumMatching


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]

    numMaximumMatching = solve(seq)

    solution = load_results(solution_path)
    correct = verify(numMaximumMatching, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_mmch_1.txt'

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]

    # min_dist = 1
    # add_adjacencies = False
    # g = bondingGraph(seq, RNA_COMPLEMENT, min_dist, add_adjacencies)
    # neato = drawBondingGraph(g)
    # neato.view()

    numMaximumMatching = solve(seq)

    print(numMaximumMatching)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(numMaximumMatching), 'w')

    # correct = solve_and_check(path)
    # print(correct)
