import networkx as nx
from BioInfoToolkit.Phylogeny.Phylo import PhyloTree
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/nkew/

In a weighted tree, each edge is assigned a (usually positive) number, called its weight. The distance between two nodes in a weighted tree becomes the sum of the weights along the unique path connecting the nodes.

To generalize Newick format to the case of a weighted tree T , during our repeated "key step," if leaves v1,v2,…,vn are neighbors in T, and all these leaves are incident to u, then we replace u with (v1:d1,v2:d2,…,vn:dn)u, where di is now the weight on the edge {vi,u}.

    Given: A collection of n weighted trees (n≤40) in Newick format, with each tree containing at most 200 nodes; each tree Tk is followed by a pair of nodes xk and yk in Tk.

    Return: A collection of n numbers, for which the kth number represents the distance between xk and yk in Tk.
"""


def verify(result: list[int], solution: list[int]) -> bool:
    correct = result == solution
    return correct


def solve(data: list[tuple[str, list[str]]]) -> list[int]:
    dists: list[int] = []
    for newick, nodePair in data:
        tree = PhyloTree(newick)
        nodeId1 = tree.findNodeByName(nodePair[0])
        nodeId2 = tree.findNodeByName(nodePair[1])
        dist = nx.algorithms.shortest_path_length(
            tree.tree, source=nodeId1, target=nodeId2, weight='weight')
        dists.append(int(dist))
    return dists


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    dists = [int(v) for v in lines[0].split()]
    return dists


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = len(lines)
    data: list[tuple[str, list[str]]] = []
    for i in range(0, n, 3):
        newick = lines[i]
        pair = lines[i+1].split(' ')
        data.append((newick, pair))

    dists = solve(data)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(dists, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_nkew_1.txt'

    lines = readTextFile(path)
    n = len(lines)

    data: list[tuple[str, list[str]]] = []
    for i in range(0, n, 3):
        newick = lines[i]
        pair = lines[i+1].split(' ')
        data.append((newick, pair))

    dists = solve(data)

    out = ' '.join(str(dist) for dist in dists)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

