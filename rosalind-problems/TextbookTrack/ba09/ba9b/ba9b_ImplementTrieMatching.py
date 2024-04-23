from collections import Counter
from BioInfoToolkit.Sequences.Trie import TrieNx
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9b/

Given a string Text and Trie(Patterns), we can quickly check whether any string from Patterns matches a prefix of Text. To do so, we start reading symbols from the beginning of Text and see what string these symbols “spell” as we proceed along the path downward from the root of the trie, as illustrated in the figure below. For each new symbol in Text, if we encounter this symbol along an edge leading down from the present node, then we continue along this edge; otherwise, we stop and conclude that no string in Patterns matches a prefix of Text. If we make it all the way to a leaf, then the pattern spelled out by this path matches a prefix of Text.

This algorithm is called PREFIXTRIEMATCHING.

    PREFIXTRIEMATCHING(Text, Trie)
        symbol ← first letter of Text
        v ← root of Trie
        while forever
            if v is a leaf in Trie
                return the pattern spelled by the path from the root to v
            else if there is an edge (v, w) in Trie labeled by symbol
                symbol ← next letter of Text
                v ← w
            else
                output "no matches found"
                return


PREFIXTRIEMATCHING finds whether any strings in Patterns match a prefix of Text. To find whether any strings in Patterns match a substring of Text starting at position k, we chop off the first k - 1 symbols from Text and run PREFIXTRIEMATCHING on the shortened string. As a result, to solve the Multiple Pattern Matching Problem (introduced in “Construct a Trie from a Collection of Patterns”, we simply iterate PREFIXTRIEMATCHING |Text| times, chopping the first symbol off of Text before each new iteration.

    TRIEMATCHING(Text, Trie)
        while Text is nonempty
            PREFIXTRIEMATCHING(Text, Trie)
            remove first symbol from Text


Implement TrieMatching

    Given: A string Text and a collection of strings Patterns.

    Return: All starting positions in Text where a string from Patterns appears as a substring.
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(text: str, patterns: list[str]) -> OutputT:
    trie = TrieNx()
    for pattern in patterns:
        trie.insert(pattern)

    # dot = trie.draw_dot()
    # dot.view()

    idxs = trie.trie_matching(text)

    return idxs


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    idxs = [int(v) for v in lines[0].split()]
    return idxs


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    text = lines[0]
    patterns = lines[1:]

    result = solve(text, patterns)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9b_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    text = lines[0]
    patterns = lines[1:]

    idxs = solve(text, patterns)

    out = ' '.join(str(i) for i in idxs)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

