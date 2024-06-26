import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
A binary tree is a tree in which each node has degree equal to at most 3. The binary tree will be our main tool in the construction of phylogenies.

A rooted tree is a tree in which one node (the root) is set aside to serve as the pinnacle of the tree. A standard graph theory exercise is to verify that for any two nodes of a tree, exactly one path connects the nodes. In a rooted tree, every node v will therefore have a single parent, or the unique node w such that the path from v to the root contains {v,w}. Any other node x adjacent to v is called a child of v because v must be the parent of x; note that a node may have multiple children. In other words, a rooted tree possesses an ordered hierarchy from the root down to its leaves, and as a result, we may often view a rooted tree with undirected edges as a directed graph in which each edge is oriented from parent to child. We should already be familiar with this idea; it's how the Rosalind problem tree works!

Even though a binary tree can include nodes having degree 2, an unrooted binary tree is defined more specifically: all internal nodes have degree 3. In turn, a rooted binary tree is such that only the root has degree 2 (all other internal nodes have degree 3).

    Given: A positive integer n (3≤n≤10000).

    Return: The number of internal nodes of any unrooted binary tree having n leaves.
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int) -> int:
    inod = n-2
    return inod


def load_results(path: str) -> int:
    lines = readTextFile(path)
    inod = int(lines[0])
    return inod


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    n = int(lines[0])
    result = solve(n)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_inod_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    inod = n - 2

    print(inod)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(inod), 'w')

    # correct = solve_and_check(path)
    # print(correct)
