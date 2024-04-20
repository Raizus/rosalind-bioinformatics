from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.GenomeAssembly import OverlapGraph
import os

"""
A graph whose nodes have all been labeled can be represented by an adjacency list, in which each row of the list contains the two node labels corresponding to a unique edge.

A directed graph (or digraph) is a graph containing directed edges, each of which has an orientation. That is, a directed edge is represented by an arrow instead of a line segment; the starting and ending nodes of an edge form its tail and head, respectively. The directed edge with tail v and head w is represented by (v,w) (but not by (w,v)). A directed loop is a directed edge of the form (v,v).

For a collection of strings and a positive integer k, the overlap graph for the strings is a directed graph Ok in which each string is represented by a node, and string s is connected to string t with a directed edge when there is a length k suffix of s that matches a length k prefix of t, as long as s≠t; we demand s≠t to prevent directed loops in the overlap graph (although directed cycles may be present).

    Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.

    Return: The adjacency list corresponding to O_3. You may return edges in any order.
"""

def verify(result: list[tuple[str, str]], solution: list[tuple[str, str]]) -> bool:
    for adj in result:
        if adj not in solution:
            return False

    return True


def solve(fasta_dict: dict[str, str], k: int) -> list[tuple[str, str]]:
    overlapGraph = OverlapGraph(fasta_dict, k)
    adj_list: list[tuple[str,str]] = []
    for edge in overlapGraph.graph.edges:
        n1, n2 = edge
        adj_list.append((n1,n2))
    return adj_list


def load_results(path: str) -> list[tuple[str,str]]:
    lines = readTextFile(path)
    adj_list = [(v1, v2) for line in lines for(v1, v2) in line.split() if len(line)]
    return adj_list


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(input_path)
    k = 3
    result = solve(fasta_dict, k)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_grph_1.txt'

    fasta_dict = read_FASTA(path)
    k = 3
    overlapGraph = OverlapGraph(fasta_dict, k)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for edge in overlapGraph.graph.edges:
        n1, n2 = edge
        print(n1 + ' ' + n2)
        writeTextFile(result_path, n1 + ' ' + n2, 'a')

    # dot = overlapGraph.draw_dot()
    # dot.view()
