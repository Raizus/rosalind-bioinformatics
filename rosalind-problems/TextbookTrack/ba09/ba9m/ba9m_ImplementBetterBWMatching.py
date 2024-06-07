from BioInfoToolkit.Sequences.BarrowsWheeler import BWTMatchingWithCheckpoints, CheckpointArrayBWT, getFirstOccurences
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9m/

If you implemented BWMATCHING in “Implement BWMatching”, you probably found the algorithm to be slow. The reason for its sluggishness is that updating the pointers top and bottom is time-intensive, since it requires examining every symbol in LastColumn between top and bottom at each step. To improve BWMATCHING, we introduce a function Countsymbol(i, LastColumn), which returns the number of occurrences of symbol in the first i positions of LastColumn. For example, Count"n”(10, "smnpbnnaaaaa$a”) = 3, and Count"a”(4, "smnpbnnaaaaa$a”) = 0.

The green lines from BWMATCHING can be compactly described without the First-to-Last mapping by the following two lines:

    top ← position of symbol with rank Countsymbol(top, LastColumn) + 1 in FirstColumn
    bottom ← position of symbol with rank Countsymbol(bottom + 1, LastColumn) in FirstColumn


Define FirstOccurrence(symbol) as the first position of symbol in FirstColumn. If Text = "panamabananas$", then FirstColumn is "$aaaaaabmnnnps", and the array holding all values of FirstOccurrence is [0, 1, 7, 8, 9, 11, 12]. For DNA strings of any length, the array FirstOccurrence contains only five elements.

The two lines of pseudocode from the previous step can now be rewritten as follows:

    top ← FirstOccurrence(symbol) + Countsymbol(top, LastColumn)
    bottom ← FirstOccurrence(symbol) + Countsymbol(bottom + 1, LastColumn) - 1


In the process of simplifying the green lines of pseudocode from BWMATCHING, we have also eliminated the need for both FirstColumn and LastToFirst, resulting in a more efficient algorithm called BETTERBWMATCHING.

    BETTERBWMATCHING(FirstOccurrence, LastColumn, Pattern, Count)
        top ← 0
        bottom ← |LastColumn| - 1
        while top ≤ bottom
            if Pattern is nonempty
                symbol ← last letter in Pattern
                remove last letter from Pattern
                if positions from top to bottom in LastColumn contain an occurrence of symbol
                    top ← FirstOccurrence(symbol) + Countsymbol(top, LastColumn)
                    bottom ← FirstOccurrence(symbol) + Countsymbol(bottom + 1, LastColumn) - 1
                else
                    return 0
            else
                return bottom - top + 1

                
Implement BetterBWMatching
    Given: A string BWT(Text), followed by a collection of strings Patterns.

    Return: A list of integers, where the i-th integer corresponds to the number of substring matches of the i-th member of Patterns in Text.
"""

OutputT = list[int]

def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(bwt: str, patterns: list[str]) -> OutputT:
    matches_count: list[int] = []
    c = 10
    checkpoint_array = CheckpointArrayBWT(bwt, c)
    first_occurence = getFirstOccurences(sorted(char for char in bwt))

    for pattern in patterns:
        match = BWTMatchingWithCheckpoints(
            bwt, pattern, checkpoint_array, first_occurence)
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
    path = f'{cwd}/rosalind_ba9m_1.txt'

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
