from BioInfoToolkit.Sequences.SequenceUtils import distanceBetweenPatternAndStrings
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba2h/

The first potential issue with implementing MedianString from “Find a Median String” is writing a function to compute d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai), the sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}. This task is achieved by the following pseudocode.

    DistanceBetweenPatternAndStrings(Pattern, Dna)
        k ← |Pattern|
        distance ← 0
        for each string Text in Dna
            HammingDistance ← ∞
            for each k-mer Pattern' in Text
                if HammingDistance > HammingDistance(Pattern, Pattern')
                    HammingDistance ← HammingDistance(Pattern, Pattern')
            distance ← distance + HammingDistance
        return distance

Compute DistanceBetweenPatternAndStrings
    Find the distance between a pattern and a set of strings.

        Given: A DNA string Pattern and a collection of DNA strings Dna.

        Return: DistanceBetweenPatternAndStrings(Pattern, Dna).

"""

OutputT = int


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(sequences: list[str], pattern: str) -> OutputT:
    dist = distanceBetweenPatternAndStrings(sequences, pattern)
    return int(dist)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    dist = int(lines[0])
    return dist


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    pattern = lines[0]
    strings = lines[1].split()

    result = solve(strings, pattern)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba2h_1.txt'

    lines = readTextFile(path)
    pattern = lines[0]
    strings = lines[1].split()

    dist = solve(strings, pattern)

    out = str(dist)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
