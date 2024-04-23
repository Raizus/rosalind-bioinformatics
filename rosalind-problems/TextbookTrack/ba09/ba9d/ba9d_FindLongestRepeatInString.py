from BioInfoToolkit.Sequences.SuffixTree import SuffixTree
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9d/

Although the suffix tree decreases memory requirements from O(|Text|2) to O(|Text|), on average it still requires about 20 times as much memory as Text. In the case of a 3 GB human genome, 60 GB of RAM still presents a memory challenge for most machines. This reveals a dark secret of big-O notation, which is that it ignores constant factors. For long strings such as the human genome, we will need to pay attention to this constant factor, since the expression O(|Text|) applies to both an algorithm with 2 · |Text| memory and an algorithm with 1000 · |Text| memory.

Yet before seeing how we can further reduce the memory needed for multiple pattern matching, we ask you to solve three problems showing how suffix trees can be applied to other pattern matching challenges. The first such problem is the Longest Repeat Problem.

Longest Repeat Problem
    Find the longest repeat in a string.

        Given: A string Text.

        Return: A longest substring of Text that appears in Text more than once. (Multiple solutions may exist, in which case you may return any one.)
"""

OutputT = str


def verify(result: OutputT, solution: OutputT, string: str) -> bool:
    correct = len(result) == len(solution) and string.find(result) != -1
    return correct


def solve(string: str) -> OutputT:
    sTree = SuffixTree(string+"$")
    # dot = sTree.drawDot()
    # dot.view()

    substrs = sTree.findLongestRepeatingSubstring(2)
    substr = substrs[0]
    return substr


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    substr = lines[0]
    return substr


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    string = lines[0]

    result = solve(string)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution, string)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9d_1.txt'

    lines = readTextFile(path)
    string = lines[0]

    substr = solve(string)

    out = substr
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

