import os
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from collections import Counter
from typing import Iterable
from BioInfoToolkit.Sequences.SequenceUtils import bondingGraph, countNoncrossingPerfectMatches
from BioInfoToolkit.Sequences.structures import RNA_COMPLEMENT
import networkx as nx

from BioInfoToolkit.IO import read_FASTA


"""
https://rosalind.info/problems/cat/

A matching in a graph is noncrossing if none of its edges cross each other. If we assume that the n nodes of this graph are arranged around a circle, and if we label these nodes with positive integers between 1 and n, then a matching is noncrossing as long as there are not edges {i,j} and {k,l} such that i<k<j<l.

A noncrossing matching of basepair edges in the bonding graph corresponding to an RNA string will correspond to a possible secondary structure of the underlying RNA strand that lacks pseudoknots, as shown in Figure 3.

In this problem, we will consider counting noncrossing perfect matchings of basepair edges. As a motivating example of how to count noncrossing perfect matchings, let cn denote the number of noncrossing perfect matchings in the complete graph K2n. After setting c0=1, we can see that c1 should equal 1 as well. As for the case of a general n, say that the nodes of K2n are labeled with the positive integers from 1 to 2n. We can join node 1 to any of the remaining 2n-1 nodes; yet once we have chosen this node (say m), we cannot add another edge to the matching that crosses the edge {1,m}. As a result, we must match all the edges on one side of {1,m} to each other. This requirement forces m to be even, so that we can write m=2k for some positive integer k.

There are 2k-2 nodes on one side of {1,m} and 2n-2k nodes on the other side of {1,m}, so that in turn there will be ck-1⋅cn-k different ways of forming a perfect matching on the remaining nodes of K2n. If we let m vary over all possible n-1 choices of even numbers between 1 and 2n, then we obtain the recurrence relation cn=∑nk=1ck-1⋅cn-k. The resulting numbers cn counting noncrossing perfect matchings in K2n are called the Catalan numbers, and they appear in a huge number of other settings. See Figure 4 for an illustration counting the first four Catalan numbers.

    Given: An RNA string s having the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'. The length of the string is at most 300 bp.

    Return: The total number of noncrossing perfect matchings of basepair edges in the bonding graph of s, modulo 1,000,000.
"""


def getIndexes(g: nx.Graph, target_nt: str) -> list[int]:
    idxs = [node for node, nattr in g.nodes.items() if nattr['nt']
            == target_nt]
    return idxs


def get_nts(g: nx.Graph, idxs: Iterable[int]) -> str:
    res = ''.join([g.nodes[node]['nt'] for node in idxs])
    return res


def canBeMatchedByCount(seq: str) -> bool:
    counts = Counter(seq)
    return counts.get('A', 0) == counts.get('U', 0) and counts.get('G', 0) == counts.get('C', 0)


def countNoncrossingPerfectMatchesGraph(g: nx.Graph, mod: int = 1000000, cache: dict[str, int] = dict()) -> int:
    n = g.number_of_nodes()

    seq = get_nts(g, sorted(g.nodes.keys()))
    if seq in cache:
        return cache[seq]

    catalanNum = 0
    if n == 0:
        catalanNum = 1
        cache[seq] = catalanNum
        return catalanNum
    if n % 2 != 0:
        cache[seq] = catalanNum
        return catalanNum
    if not canBeMatchedByCount(seq):
        cache[seq] = catalanNum
        return catalanNum
     
    orderedNodes: list[int] = list(sorted(g.nodes.keys()))
    node1 = orderedNodes[0]
    adjNodes = g.adj[node1].keys()
    lastNode = orderedNodes[-1]
    # note node2 > node1 always
    
    for node2 in adjNodes:
        iidxs = list(range(node1+1, node2))
        oidxs = list(range(node2+1, lastNode+1))

        subgraph1 = nx.subgraph(g, iidxs)
        c1 = countNoncrossingPerfectMatchesGraph(subgraph1, mod)
        subgraph2 = nx.subgraph(g, oidxs)
        c2 = countNoncrossingPerfectMatchesGraph(subgraph2, mod)

        catalanNum = (catalanNum + c1 * c2) % mod

    cache[seq] = catalanNum

    return catalanNum


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(seq: str) -> int:
    catalanNum = countNoncrossingPerfectMatches(seq, RNA_COMPLEMENT, mod)
    return catalanNum


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]

    mod = 1000000
    catalanNum = countNoncrossingPerfectMatches(
        seq, RNA_COMPLEMENT, min_dist, mod)

    solution = load_results(solution_path)
    correct = verify(catalanNum, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_cat_1.txt'

    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]

    mod = 1000000
    min_dist = 1
    add_adjacencies = False
    g = bondingGraph(seq, RNA_COMPLEMENT, min_dist, add_adjacencies)

    # slow
    # catalanNum = countNoncrossingPerfectMatchesGraph(g, mod)

    # fast
    catalanNum = countNoncrossingPerfectMatches(seq, RNA_COMPLEMENT, min_dist, mod)

    print(catalanNum)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(catalanNum), 'w')

    correct = solve_and_check(path)
    print(correct)
