from collections import Counter
from BioInfoToolkit.Sequences.SuffixTree import SuffixTree
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/suff/

Given a string s having length n, recall that its suffix tree T(s) is defined by the following properties:

    - T(s) is a rooted tree having exactly n leaves.
    - Every edge of T(s) is labeled with a substring of s*, where s* is the string formed by adding a placeholder symbol $ to the end of s.
    - Every internal node of T(s) other than the root has at least two children; i.e., it has degree at least 3.
    - The substring labels for the edges leading down from a node to its children must begin with different symbols.
    - By concatenating the substrings along edges, each path from the root to a leaf corresponds to a unique suffix of s*.

Figure 1 contains an example of a suffix tree.

    Given: A DNA string s of length at most 1kbp.

    Return: The substrings of s* encoding the edges of the suffix tree for s. You may list these substrings in any order.
"""


def verify(result: list[str], solution: list[str]) -> bool:
    res_counter = Counter(result)
    sol_counter = Counter(solution)

    correct = res_counter == sol_counter
    return correct


def solve(seq: str) -> list[str]:
    sTree = SuffixTree(seq)
    words: list[str] = []
    for _, nbrs in sTree.graph.adj.items():
        for _, edge in nbrs.items():
            start, end = edge['start'], edge['end']
            label = sTree.seq[start:end]
            words.append(label)
    return words



def load_results(path: str) -> list[str]:
    lines = readTextFile(path)
    words = [line for line in lines if len(line) and not line.isspace()]
    return words


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq = lines[0]

    result = solve(seq)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_suff_1.txt'

    lines = readTextFile(path)
    seq = lines[0]

    sTree = SuffixTree(seq)
    # dot = sTree.drawDot()
    # dot.view()

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for n1, nbrs in sTree.graph.adj.items():
        for n2, edge in nbrs.items():
            start, end = edge['start'], edge['end']
            label = sTree.seq[start:end]
            print(label)
            writeTextFile(result_path, label, 'a')

    # correct = solve_and_check(path)
    # print(correct)
