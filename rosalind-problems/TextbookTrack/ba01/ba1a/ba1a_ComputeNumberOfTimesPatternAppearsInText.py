from BioInfoToolkit.Sequences.StringUtils import count_kmer
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba1a/

This is the first problem in a collection of "code challenges" to accompany Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.

A k-mer is a string of length k. We define Count(Text, Pattern) as the number of times that a k-mer Pattern appears as a substring of Text. For example,

    Count(ACAACTATGCATACTATCGGGAACTATCCT,ACTAT)=3.

We note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2) since we should account for overlapping occurrences of Pattern in Text.

To compute Count(Text, Pattern), our plan is to “slide a window” down Text, checking whether each k-mer substring of Text matches Pattern. We will therefore refer to the k-mer starting at position i of Text as Text(i, k). Throughout this book, we will often use 0-based indexing, meaning that we count starting at 0 instead of 1. In this case, Text begins at position 0 and ends at position |Text| - 1 (|Text| denotes the number of symbols in Text). For example, if Text = GACCATACTG, then Text(4, 3) = ATA. Note that the last k-mer of Text begins at position |Text| - k, e.g., the last 3-mer of GACCATACTG starts at position 10 - 3 = 7. This discussion results in the following pseudocode for computing Count(Text, Pattern).

    PatternCount(Text, Pattern)
        count ← 0
        for i ← 0 to |Text| - |Pattern|
            if Text(i, |Pattern|) = Pattern
                count ← count + 1
        return count

    Given: {DNA strings}} Text and Pattern.

    Return: Count(Text, Pattern).
"""


def verify(result: int, solution: int) -> bool:    
    correct = result == solution
    return correct


def solve(seq: str, kmer: str) -> int:
    count = count_kmer(string, pattern)
    return count


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    string = lines[0]
    pattern = lines[1]
    result = solve(string, pattern)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba1a_1.txt'

    lines = readTextFile(path)
    string = lines[0]
    pattern = lines[1]

    count = count_kmer(string, pattern)

    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    correct = solve_and_check(path)
    print(correct)
