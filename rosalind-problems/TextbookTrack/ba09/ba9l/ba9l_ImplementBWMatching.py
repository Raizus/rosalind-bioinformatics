from BioInfoToolkit.Sequences.BarrowsWheeler import BWTMatchingWithCheckpoints
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9l/

We are now ready to describe BWMATCHING, an algorithm that counts the total number of matches of Pattern in Text, where the only information that we are given is FirstColumn and LastColumn = BWT(Text) in addition to the Last-to-First mapping (from “Generate the Last-to-First Mapping of a String”). The pointers top and bottom are updated by the green lines in the following pseudocode.

    BWMATCHING(FirstColumn, LastColumn, Pattern, LastToFirst)
        top ← 0
        bottom ← |LastColumn| - 1
        while top ≤ bottom
            if Pattern is nonempty
                symbol ← last letter in Pattern
                remove last letter from Pattern
                if positions from top to bottom in LastColumn contain an occurrence of symbol
                    topIndex ← first position of symbol among positions from top to bottom in LastColumn
                    bottomIndex ← last position of symbol among positions from top to bottom in LastColumn
                    top ← LastToFirst(topIndex)
                    bottom ← LastToFirst(bottomIndex)
                else
                    return 0
            else
                return bottom - top + 1

Implement BWMatching

    Given: A string BWT(Text), followed by a collection of strings Patterns.

    Return: A list of integers, where the i-th integer corresponds to the number of substring matches of the i-th member of Patterns in Text.
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(bwt: str, patterns: list[str]) -> OutputT:
    matches_count: list[int] = []
    for pattern in patterns:
        match = BWTMatchingWithCheckpoints(bwt, pattern)
        if match is not None:
            top, bottom = match
            n = bottom - top + 1
            matches_count.append(n)
        else:
            matches_count.append(0)
    return matches_count


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    matches_count = [int(v) for v in lines[0].split()]
    return matches_count


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    bwt = lines[0]
    patterns = lines[1].split()

    result = solve(bwt, patterns)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9l_1.txt'

    lines = readTextFile(path)
    bwt = lines[0]
    patterns = lines[1].split()

    matches_count = solve(bwt, patterns)

    out = ' '.join(str(v) for v in matches_count)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
