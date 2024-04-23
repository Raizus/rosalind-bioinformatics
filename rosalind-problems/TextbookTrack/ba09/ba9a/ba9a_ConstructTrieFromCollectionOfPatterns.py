from BioInfoToolkit.Sequences.Trie import TrieNx
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import re


"""
https://rosalind.info/problems/ba9a/

Reads will form a collection of strings Patterns that we wish to match against a reference genome Text. For each string in Patterns, we will first find all its exact matches as a substring of Text (or conclude that it does not appear in Text). When hunting for the cause of a genetic disorder, we can immediately eliminate from consideration areas of the reference genome where exact matches occur. We will later generalize this problem to find approximate matches, where single nucleotide substitutions in reads separate the individual from the reference genome (or represent errors in reads).

Multiple Pattern Matching Problem: Find all occurrences of a collection of patterns in a text.
    Input: A string Text and a collection Patterns containing (shorter) strings.
    Output: All starting positions in Text where a string from Patterns appears as a substring.

To solve this problem, we will consolidate Patterns into a directed tree called a trie (pronounced “try”), which is written Trie(Patterns) and has the following properties.

    - The trie has a single root node with indegree 0, denoted root.
    - Each edge of Trie(Patterns) is labeled with a letter of the alphabet.
    - Edges leading out of a given node have distinct labels.
    - Every string in Patterns is spelled out by concatenating the letters along some path from the root downward.
    - Every path from the root to a leaf, or node with outdegree 0, spells a string from Patterns.

The most obvious way to construct Trie(Patterns) is by iteratively adding each string from Patterns to the growing trie, as implemented by the following algorithm.

    TRIECONSTRUCTION(Patterns)
        Trie ← a graph consisting of a single node root
        for each string Pattern in Patterns
            currentNode ← root
            for i ← 1 to |Pattern|
                if there is an outgoing edge from currentNode with label currentSymbol
                    currentNode ← ending node of this edge
                else
                    add a new node newNode to Trie
                    add a new edge from currentNode to newNode with label currentSymbol
                    currentNode ← newNode
        return Trie


Trie Construction Problem
    Construct a trie on a collection of patterns.

        Given: A collection of strings Patterns.

        Return: The adjacency list corresponding to Trie(Patterns), in the following format. If Trie(Patterns) has n nodes, first label the root with 1 and then label the remaining nodes with the integers 2 through n in any order you like. Each edge of the adjacency list of Trie(Patterns) will be encoded by a triple: the first two members of the triple must be the integers labeling the initial and terminal nodes of the edge, respectively; the third member of the triple must be the symbol labeling the edge.
"""

OutputT = list[tuple[int,int,str]]


def verify(result: OutputT, solution: OutputT) -> bool:
    # TODO: check if correct
    return False


def solve(patterns: list[str]) -> OutputT:
    trie = TrieNx()
    for pattern in patterns:
        trie.insert(pattern)

    # dot = trie.draw_dot()
    # dot.view()

    adj_list: list[tuple[int, int, str]] = []

    for edge in trie.graph.edges:
        adj_list.append((edge[0], edge[1], trie.graph.edges[edge]['char']))

    return adj_list


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    adj_list: OutputT = []
    for line in lines:
        if not len(line) or line.isspace():
            continue

        matches = re.match(r'(\d+)->(\d+):(\w)', line)
        if matches is None:
            raise Exception('Line {line} does not have the correct format')

        n1, n2, label = int(matches[1]), int(matches[2]), matches[3]
        adj_list.append((n1, n2, label))

    return adj_list


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seqs = [line for line in lines if len(line) and not line.isspace()]

    result = solve(seqs)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9a_1.txt'

    lines = readTextFile(path)
    seqs = [line for line in lines if len(line) and not line.isspace()]

    adj_list = solve(seqs)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for n1, n2, char in adj_list:
        out = f"{n1}->{n2}:{char}"
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
