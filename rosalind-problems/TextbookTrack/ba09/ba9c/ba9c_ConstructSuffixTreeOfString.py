from BioInfoToolkit.Sequences.SuffixTree import SuffixTree
from collections import Counter
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9c/

Storing Trie(Patterns), which was introduced in “Implement TrieMatching”, requires a great deal of memory. So let's process Text into a data structure instead. Our goal is to compare each string in Patterns against Text without needing to to traverse Text from beginning to end. In more familiar terms, instead of packing Patterns onto a bus and riding the long distance down Text, our new data structure will be able to “teleport” each string in Patterns directly to its occurrences in Text.

A suffix trie, denoted SuffixTrie(Text), is the trie formed from all suffixes of Text. From now on, we append the dollar-sign ("$") to Text in order to mark the end of Text. We will also label each leaf of the resulting trie by the starting position of the suffix whose path through the trie ends at this leaf (using 0-based indexing). This way, when we arrive at a leaf, we will immediately know where this suffix came from in Text.

However, the runtime and memory required to construct SuffixTrie(Text) are both equal to the combined length of all suffixes in Text. There are |Text| suffixes of Text, ranging in length from 1 to |Text| and having total length |Text|·(|Text|+1)/2, which is O(|Text|2). Thus, we need to reduce both the construction time and memory requirements of suffix tries to make them practical.

Let's not give up hope on suffix tries. We can reduce the number of edges in SuffixTrie(Text) by combining the edges on any non-branching path into a single edge. We then label this edge with the concatenation of symbols on the consolidated edges, as shown in the figure below. The resulting data structure is called a suffix tree, written SuffixTree(Text).

To match a single Pattern to Text, we thread Pattern into SuffixTree(Text) by the same process used for a suffix trie. Similarly to the suffix trie, we can use the leaf labels to find starting positions of successfully matched patterns.

Suffix trees save memory because they do not need to store concatenated edge labels from each non-branching path. For example, a suffix tree does not need ten bytes to store the edge labeled "mabananas$" in SuffixTree("panamabananas$") (reproduced below); instead, it suffices to store a pointer to position 4 of "panamabananas$", as well as the length of "mabananas$". Furthermore, suffix trees can be constructed in linear time, without having to first construct the suffix trie! We will not ask you to implement this fast suffix tree construction algorithm because it is quite complex.

Suffix Tree Construction Problem
    Construct the suffix tree of a string.

        Given: A string Text.

        Return: The strings labeling the edges of SuffixTree(Text). (You may return these strings in any order.)
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(string: str) -> OutputT:
    sTree = SuffixTree(string)
    words: list[str] = []
    for _, nbrs in sTree.graph.adj.items():
        for _, edge in nbrs.items():
            start, end = edge['start'], edge['end']
            label = sTree.seq[start:end]
            words.append(label)
    return words


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    strings = [line for line in lines if len(line) and not line.isspace()]
    return strings


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    string = lines[0]

    result = solve(string)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9c_1.txt'

    lines = readTextFile(path)
    string = lines[0]

    words = solve(string)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for word in words:
        print(word)
        writeTextFile(result_path, word, 'a')

    correct = solve_and_check(path)
    print(correct)


