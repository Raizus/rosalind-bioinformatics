
import networkx as nx
from BioInfoToolkit.Phylogeny.Phylo import PhyloTree
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/nwck/

Newick format is a way of representing trees even more concisely than using an adjacency list, especially when dealing with trees whose internal nodes have not been labeled.

First, consider the case of a rooted tree T. A collection of leaves v1,v2,…,vn of T are neighbors if they are all adjacent to some internal node u. Newick format for T is obtained by iterating the following key step: delete all the edges {vi,u} from T and label u with (v1,v2,…,vn)u. This process is repeated all the way to the root, at which point a semicolon signals the end of the tree.

A number of variations of Newick format exist. First, if a node is not labeled in T, then we simply leave blank the space occupied by the node. In the key step, we can write (v1,v2,…,vn) in place of (v1,v2,…,vn)u if the vi are labeled; if none of the nodes are labeled, we can write (,,…,).

A second variation of Newick format occurs when Tis unrooted, in which case we simply select any internal node to serve as the root of T. A particularly peculiar case of Newick format arises when we choose a leaf to serve as the root.

Note that there will be a large number of different ways to represent T in Newick format; see Figure 1.

    Given: A collection of n trees (n≤40) in Newick format, with each tree containing at most 200 nodes; each tree Tk is followed by a pair of nodes xk and yk in Tk.

    Return: A collection of n positive integers, for which the kth integer represents the distance between xk and yk in Tk.
"""


def verify(result: list[int], solution: list[int]) -> bool:
    correct = result == solution
    return correct


def solve(data: list[tuple[str, list[str]]]) -> list[int]:
    dists = []
    for newick, nodePair in data:
        tree = PhyloTree(newick)
        nodeId1 = tree.findNodeByName(nodePair[0])
        nodeId2 = tree.findNodeByName(nodePair[1])
        dist = nx.algorithms.shortest_path_length(
            tree.tree, source=nodeId1, target=nodeId2)
        dists.append(dist)
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
    path = f'{cwd}/rosalind_nwck_1.txt'

    lines = readTextFile(path)
    n = len(lines)
    data: list[tuple[str, list[str]]] = []
    for i in range(0, n, 3):
        newick = lines[i]
        pair = lines[i+1].split(' ')
        data.append((newick, pair))

    dists = solve(data)

    # tree = parseNewick(test)
    # dot = draw_dot(tree)
    # dot.view()

    out = ' '.join(str(dist) for dist in dists)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
