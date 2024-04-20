from collections import defaultdict
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
An undirected graph is connected if there is a path connecting any two nodes. A tree is a connected (undirected) graph containing no cycles; this definition forces the tree to have a branching structure organized around a central core of nodes, just like its living counterpart. See Figure 2.

We have already grown familiar with trees in “Mendel's First Law”, where we introduced the probability tree diagram to visualize the outcomes of a random variable.

In the creation of a phylogeny, taxa are encoded by the tree's leaves, or nodes having degree 1. A node of a tree having degree larger than 1 is called an internal node.

    Given: A positive integer n (n≤1000) and an adjacency list corresponding to a graph on n nodes that contains no cycles.

    Return: The minimum number of edges that can be added to the graph to produce a tree.
"""


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, edges: list[tuple[int, ...]]) -> int:
    graph: defaultdict[int, set[int]] = defaultdict(set)
    for n1, n2 in edges:
        graph[n1].add(n2)
        graph[n2].add(n1)

        # To count the minimum number of edges that can be added to the graph to produce a tree is equal to the number of islands - 1
    island_count = 0

    # choose starting node
    visited: set[int] = set()
    tracking_nodes: set[int] = set()
    node0 = next(iter(graph))
    tracking_nodes.add(node0)

    while len(graph):
        node0 = next(iter(graph))
        tracking_nodes.add(node0)
        while len(tracking_nodes):
            temp: set[int] = set()
            for node in tracking_nodes:
                if node not in graph:
                    continue
                visited.add(node)
                adj = graph.pop(node)
                temp.update(adj)

            temp.difference_update(visited)
            tracking_nodes.clear()
            tracking_nodes.update(temp)
        island_count += 1

    num_needed_edge = island_count - 1 + num_nodes-len(visited)
    return num_needed_edge


def load_results(path: str) -> int:
    lines = readTextFile(path)
    n_edges = int(lines[0])
    return n_edges


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(path)
    num_nodes = int(lines[0])
    edges: list[tuple[int, ...]] = [
        tuple(int(n) for n in line.split()) for line in lines[1:]]

    result = solve(num_nodes, edges)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_tree_1.txt'

    lines = readTextFile(path)
    num_nodes = int(lines[0])
    edges: list[tuple[int, ...]] = [tuple(int(n) for n in line.split()) for line in lines[1:]]

    num_needed_edges = solve(num_nodes, edges)

    print(num_needed_edges)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(num_needed_edges), 'w')
