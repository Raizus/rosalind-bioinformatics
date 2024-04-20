from BioInfoToolkit.Sequences.SuffixTree import SuffixTree
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/mrep/

A maximal repeat of a string s is a repeated substring t of s having two occurrences t1 and t2 such that t1 and t2 cannot be extended by one symbol in either direction in s and still agree.

For example, "AG" is a maximal repeat in "TAGTTAGCGAGA" because even though the first two occurrences of "AG" can be extended left into "TAG", the first and third occurrences differ on both sides of the repeat; thus, we conclude that "AG" is a maximal repeat. Note that "TAG" is also a maximal repeat of "TAGTTAGCGAGA", since its only two occurrences do not still match if we extend them in either direction.

    Given: A DNA string s of length at most 1 kbp.

    Return: A list containing all maximal repeats of s having length at least 20.
"""

def verify(result: list[str], solution: list[str]) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(seq: str, k: int) -> list[str]:
    sTree = SuffixTree(seq+'$')
    maximal_repeats = sTree.findMaximalRepeats(k)
    return maximal_repeats


def load_results(path: str) -> list[str]:
    lines = readTextFile(path)
    strings = [line for line in lines if len(line) and not line.isspace()]
    return strings


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq = lines[0]
    k = 20

    maximal_repeats = solve(seq, k)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(maximal_repeats, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_mrep_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    k = 20

    sTree = SuffixTree(seq+'$')
    maximalRepeats = sTree.findMaximalRepeats(k)
    # dot = sTree.drawDot()
    # dot.view()

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for a in maximalRepeats:
        print(a)
        writeTextFile(result_path, a, 'a')

    correct = solve_and_check(path)
    print(correct)
