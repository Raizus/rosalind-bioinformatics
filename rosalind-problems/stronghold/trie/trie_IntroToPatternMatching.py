
from BioInfoToolkit.Sequences.Trie import Trie, TrieNode
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.IO import readTextFile, writeTextFile
import re


"""
https://rosalind.info/problems/trie/

Given a collection of strings, their trie (often pronounced "try" to avoid ambiguity with the general term tree) is a rooted tree formed as follows. For every unique first symbol in the strings, an edge is formed connecting the root to a new vertex. This symbol is then used to label the edge.

We may then iterate the process by moving down one level as follows. Say that an edge connecting the root to a node v is labeled with 'A'; then we delete the first symbol from every string in the collection beginning with 'A' and then treat v as our root. We apply this process to all nodes that are adjacent to the root, and then we move down another level and continue. See Figure 1 for an example of a trie.

As a result of this method of construction, the symbols along the edges of any path in the trie from the root to a leaf will spell out a unique string from the collection, as long as no string is a prefix of another in the collection (this would cause the first string to be encoded as a path terminating at an internal node).

    Given: A list of at most 100 DNA strings of length at most 100 bp, none of which is a prefix of another.

    Return: The adjacency list corresponding to the trie T for these patterns, in the following format. If T has n nodes, first label the root with 1 and then label the remaining nodes with the integers 2 through n in any order you like. Each edge of the adjacency list of T will be encoded by a triple containing the integer representing the edge's parent node, followed by the integer representing the edge's child node, and finally the symbol labeling the edge.
"""


def save_result(result_path: str, node: TrieNode):
    def dfs(node: TrieNode):
        for label, childNode in node.children.items():
            res = f"{node.id+1} {childNode.id+1} {label}"
            print(res)
            writeTextFile(result_path, res, "a")
            dfs(childNode)

    dfs(node)
    

def verify(result: list[tuple[int, int, str]], seqs: list[str]) -> bool:

    return False


def solve(seqs: list[str]) -> list[tuple[int, int, str]]:
    trie = Trie()
    for seq in seqs:
        trie.insert(seq)

    adj_list: list[tuple[int, int, str]] = []

    def dfs(node: TrieNode):
        for label, childNode in node.children.items():
            adj_list.append((node.id, childNode.id, str(label)))

    dfs(trie.root)
    return adj_list


def load_results(path: str) -> list[tuple[int,int,str]]:
    lines = readTextFile(path)
    adj_list: list[tuple[int, int, str]] = []
    for line in lines:
        if not len(line) or line.isspace():
            continue

        matches = re.match(r'(\d+)\s+(\d+)\s+(\w)', line)
        if matches is None:
            raise Exception('Line {line} does not have the correct format')

        n1, n2, label = int(matches[1]), int(matches[2]), matches[3]
        adj_list.append((n1,n2,label))

    return adj_list


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seqs = [line for line in lines if len(line) and not line.isspace()]
    adj_list = solve(seqs)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(adj_list, seqs)
    return correct

# TODO finish this

if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_trie_0.txt'

    lines = readTextFile(path)
    seqs = [line for line in lines if len(line) and not line.isspace()]
    
    trie = Trie()
    for seq in seqs:
        trie.insert(seq)

    dot = trie.draw_dot()
    # dot.view()
    dot.render('trie', directory=f'{cwd}/')

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    save_result(result_path, trie.root)

    # correct = solve_and_check(path)
    # print(correct)
