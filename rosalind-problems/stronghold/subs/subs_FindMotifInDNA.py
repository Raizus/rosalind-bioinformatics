from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.BioSequence import BioSequence
import os

"""
Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s (as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i of s is denoted by s[i].

A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s.

    Given: Two DNA strings s and t (each of length at most 1 kbp).

    Return: All locations of t as a substring of s.
"""


def verify(result: list[int], solution: list[int]) -> bool:
    correct = all(r == s for r,s in zip(result, solution))
    return correct


def solve(seq: str, pattern: str) -> list[int]:
    bio_seq = BioSequence(seq)
    result = bio_seq.find_substring(pattern)
    return result


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    idxs = [int(v) for v in lines[0].split()]
    return idxs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    seq = lines[0]
    pattern = lines[1]

    idxs = solve(seq, pattern)
    idxs = [idx+1 for idx in idxs]

    solution = load_results(solution_path)

    correct = verify(idxs, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_subs_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    pattern = lines[1]

    idxs = solve(seq, pattern)
    out = ' '.join(str(idx+1) for idx in idxs)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')
