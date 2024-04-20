
import os
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.SequenceUtils import bondingGraph, countNoncrossingMatches
from BioInfoToolkit.Sequences.structures import RNA_COMPLEMENT

from BioInfoToolkit.IO import read_FASTA

"""
Similarly to our definition of the Catalan numbers, the n-th Motzkin number mn counts the number of ways to form a (not necessarily perfect) noncrossing matching in the complete graph Kn containing n nodes. For example, Figure 1 demonstrates that m5=21. Note in this figure that technically, the "trivial" matching that contains no edges at all is considered to be a matching, because it satisfies the defining condition that no two edges are incident to the same node.

How should we compute the Motzkin numbers? As with Catalan numbers, we will take m0=m1=1. To calculate mn in general, assume that the nodes of Kn are labeled around the outside of a circle with the integers between 1 and n, and consider node 1, which may or may not be involved in a matching. If node 1 is not involved in a matching, then there are mn-1 ways of matching the remaining n-1 nodes. If node 1 is involved in a matching, then say it is matched to node k: this leaves k-2 nodes on one side of edge {1,k} and n-k nodes on the other side; as with the Catalan numbers, no edge can connect the two sides, which gives us mk-2⋅mn-k ways of matching the remaining edges. Allowing k to vary between 2 and n yields the following recurrence relation for the Motzkin numbers: mn=mn-1+∑nk=2mk-2⋅mn-k.

To count all possible secondary structures of a given RNA string that do not contain pseudoknots, we need to modify the Motzkin recurrence so that it counts only matchings of basepair edges in the bonding graph corresponding to the RNA string; see Figure 2.

    Given: An RNA string s of length at most 300 bp.

    Return: The total number of noncrossing matchings of basepair edges in the bonding graph of s, modulo 1,000,000.
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq: str, min_dist: int, mod: int) -> int:
    motzkin_num = countNoncrossingMatches(seq, RNA_COMPLEMENT, min_dist, mod)
    return motzkin_num


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]

    mod = 1000000
    min_dist = 1
    motzkin_num = solve(seq, min_dist, mod)

    solution = load_results(solution_path)
    correct = verify(motzkin_num, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_motz_1.txt'

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]
    mod = 10**6

    min_dist = 1
    add_adjacencies = False
    # g = bondingGraph(
    #     seq, RNA_COMPLEMENT, min_dist, add_adjacencies)

    motzkinNum = countNoncrossingMatches(seq, RNA_COMPLEMENT, min_dist, mod)
    print(motzkinNum)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(motzkinNum), 'w')

    correct = solve_and_check(path)
    print(correct)
