from BioInfoToolkit.Sequences.BarrowsWheeler import BWTMultiplePatternMatching, barrowsWheelerTransform
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba9n/

We hope that you have noticed what is still a limitation of BETTERBWMATCHING from “Implement BetterBWMatching” — even though this algorithm counts the number of occurrences of Pattern in Text, it does not tell us where these occurrences are located in Text! To locate pattern matches identified by BETTERBWMATCHING, we can once again use the suffix array.

However, recall that our original motivation for using the Burrows-Wheeler transform was to reduce the amount of memory used by the suffix array for pattern matching. If we add the suffix array to Burrows-Wheeler-based pattern matching, then we are right back where we started!

The memory-saving device that we will employ is inelegant but useful. We will build a partial suffix array of Text, denoted SuffixArrayK(Text), which only contains values that are multiples of some positive integer K (see “Construct the Partial Suffix Array of a String”). In real applications, partial suffix arrays are often constructed for K = 100, thus reducing memory usage by a factor of 100 compared to a full suffix array.

We should also discuss how to improve BETTERBWMATCHING (reproduced below) by resolving the trade-off between precomputing the values of Countsymbol(i, LastColumn) (requiring substantial memory) and computing these values as we go (requiring substantial runtime).

The balance that we strike is similar to the one used for the partial suffix array. Rather than storing Countsymbol(i, LastColumn) for all positions i, we will only store the Count arrays when i is divisible by C, where C is a constant; these arrays are called checkpoint arrays. When C is large (C is typically equal to 100 in practice) and the alphabet is small (e.g., 4 nucleotides), checkpoint arrays require only a fraction of the memory used by BWT(Text).

What about runtime? Using checkpoint arrays, we can compute the top and bottom pointers in a constant number of steps (i.e., fewer than C). Since each Pattern requires at most |Pattern| pointer updates, the modified BETTERBWMATCHING algorithm now requires O(|Patterns|) runtime, which is the same as using a trie or suffix array.

Furthermore, we now only have to store the following data in memory: BWT(Text), FirstOccurrence, the partial suffix array, and the checkpoint arrays. Storing this data requires memory approximately equal to 1.5 · |Text|. Thus, we have finally knocked down the memory required for solving the Multiple Pattern Matching Problem for millions of sequencing reads into a workable range, and we can at last solve a full version of the Multiple Pattern Matching Problem.

Multiple Pattern Matching Problem
    Find all occurrences of a collection of patterns in a text.

        Given: A string Text and a collection of strings Patterns.

        Return: All starting positions in Text where a string from Patterns appears as a substring.
"""

OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(bwt: str, patterns: list[str]) -> OutputT:
    bwt = barrowsWheelerTransform(string+'$')
    c = 10
    matches = BWTMultiplePatternMatching(bwt, patterns, c)
    return sorted(matches)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    matches_count = [int(v) for v in lines[0].split()]
    return matches_count


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    string = lines[0]
    patterns = lines[1:]

    result = solve(string, patterns)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9n_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    string = lines[0]
    patterns = lines[1:]

    matches_count = solve(string, patterns)

    out = ' '.join(str(v) for v in matches_count)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

