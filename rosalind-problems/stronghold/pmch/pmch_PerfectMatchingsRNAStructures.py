

from collections import Counter
import math
from BioInfoToolkit.Sequences.SequenceUtils import bondingGraph, drawBondingGraph
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
A matching in a graph G is a collection of edges of G for which no node belongs to more than one edge in the collection. See Figure 2 for examples of matchings. If G contains an even number of nodes (say 2n), then a matching on G is perfect if it contains n edges, which is clearly the maximum possible. An example of a graph containing a perfect matching is shown in Figure 3.

First, let Kn denote the complete graph on 2n labeled nodes, in which every node is connected to every other node with an edge, and let pn denote the total number of perfect matchings in Kn. For a given node x, there are 2n−1 ways to join x to the other nodes in the graph, after which point we must form a perfect matching on the remaining 2n−2 nodes. This reasoning provides us with the recurrence relation pn=(2n−1)⋅pn−1; using the fact that p1 is 1, this recurrence relation implies the closed equation pn=(2n−1)(2n−3)(2n−5)⋯(3)(1).

Given an RNA string s=s1…sn, a bonding graph for s is formed as follows. First, assign each symbol of s to a node, and arrange these nodes in order around a circle, connecting them with edges called adjacency edges. Second, form all possible edges {A, U} and {C, G}, called basepair edges; we will represent basepair edges with dashed edges, as illustrated by the bonding graph in Figure 4.

Note that a matching contained in the basepair edges will represent one possibility for base pairing interactions in s, as shown in Figure 5. For such a matching to exist, s must have the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.

    Given: An RNA string s of length at most 80 bp having the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.

    Return: The total possible number of perfect matchings of basepair edges in the bonding graph of s.
"""


def verify(result: int, solution: int) -> bool:
    return result == solution


def solve(seq: str) -> int:
    counts = Counter(seq)
    if counts['A'] == counts['U'] and counts['G'] == counts['C']:
        # The number of perfect matches will be equal to the all the possible ways to match all the A's to all the U's times the number of possible ways to match all the G's to all the C's
        numPmatches = math.factorial(counts['A']) * math.factorial(counts['G'])
        return numPmatches
    raise Exception("For a perfect match to exist, the seuqence must have the same number of 'A' and 'U' and the same number of 'G' and 'C'")


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]
    result = solve(seq)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_pmch_1.txt'

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]

    # min_dist = 1
    # add_adjacencies = True

    # G = bondingGraph(
    #     seq, RNA_COMPLEMENT, min_dist, add_adjacencies)

    # neato = drawBondingGraph(G)
    # neato.view()

    counts = solve(seq)
    print(counts)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(counts), 'w')

    correct = solve_and_check(path)
    print(correct)
