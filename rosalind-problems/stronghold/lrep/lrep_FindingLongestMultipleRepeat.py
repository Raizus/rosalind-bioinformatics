from BioInfoToolkit.Sequences.SuffixTree import SuffixTree
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import re

"""
https://rosalind.info/problems/lrep/

A repeated substring of a string s of length n is simply a substring that appears in more than one location of s; more specifically, a k-fold substring appears in at least k distinct locations.

The suffix tree of s, denoted T(s), is defined as follows:

    - T(s) is a rooted tree having exactly n leaves.
    - Every edge of T(s) is labeled with a substring of s*, where s* is the string formed by adding a placeholder symbol $ to the end of s.
    - Every internal node of T(s) other than the root has at least two children; i.e., it has degree at least 3.
    - The substring labels for the edges leading from a node to its children must begin with different symbols.
    - By concatenating the substrings along edges, each path from the root to a leaf corresponds to a unique suffix of s*.

See Figure 1 for an example of a suffix tree.

    Given: A DNA string s (of length at most 20 kbp) with $ appended, a positive integer k, and a list of edges defining the suffix tree of s. Each edge is represented by four components:

        1. the label of its parent node in T(s);
        2. the label of its child node in T(s);
        3. the location of the substring t of s* assigned to the edge; and
        4. the length of t.

    Return: The longest substring of s that occurs at least k times in s. (If multiple solutions exist, you may return any single solution.)
"""

def count_occurences(seq: str, pattern: str) -> int:
    matches = re.findall(f'(?=({pattern}))', seq)
    count = len(matches)
    return count

def verify(result: str, solution: str, seq: str) -> bool:
    count_res = count_occurences(seq, result)
    count_sol = count_occurences(seq, solution)
    correct = len(result) == len(solution) and count_res == count_sol
    return correct


def solve(seq: str, k: int) -> str:
    sTree = SuffixTree(seq)
    substrs = sTree.findLongestRepeatingSubstring(k)
    return substrs[0]


def load_results(path: str) -> str:
    lines = readTextFile(path)
    longset_substring = lines[0]
    return longset_substring


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq = lines[0]
    k = int(lines[1])

    result = solve(seq, k)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, seq)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_lrep_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    k = int(lines[1])

    sTree = SuffixTree(seq)
    # dot = sTree.drawDot()
    # dot.view()

    substrs = sTree.findLongestRepeatingSubstring(k)
    print(substrs[0])

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, substrs[0], 'w')

    # correct = solve_and_check(path)
    # print(correct)
