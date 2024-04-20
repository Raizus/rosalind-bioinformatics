from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/itwv/

Given three strings s, t, and u, we say that t and u can be interwoven into s if there is some substring of s made up of t and u as disjoint subsequences.

For example, the strings "ACAG " and "CCG" can be interwoven into "GACCACGGTT". However, they cannot be interwoven into "GACCACAAAAGGTT" because of the appearance of the four 'A's in the middle of the subsequences. Similarly, even though both "ACACG" is a shortest common supersequence of ACAG and CCG, it is not possible to interweave these two strings into "ACACG" because the two desired subsequences must be disjoint; see “Interleaving Two Motifs” for details on finding a shortest common supersequence of two strings.

    Given: A text DNA string s of length at most 10 kbp, followed by a collection of n (n≤10) DNA strings of length at most 10 bp acting as patterns.

    Return: An n x n matrix M for which Mj,k=1 if the jth and kth pattern strings can be interwoven into s and Mj,k=0 otherwise.
"""


def canInterweave(pat: str, s: str, t: str) -> bool:
    """Returns True if two motifs s and t can be interwoven into a substring of pattern

    Args:
        pat (str): pattern
        s (str): motif 1
        t (str): motif 2
    """
    lpat, ls, lt = len(pat), len(s), len(t)

    def test(x: int, i: int, j: int):
        if i == ls and j == lt:
            return True
        if lpat == x:
            return False

        v = False
        if i < ls and pat[x] == s[i]:
            v = test(x+1, i+1, j)
        if not v and j < lt and pat[x] == t[j]:
            v = test(x+1, i, j+1)
        return v

    for i in range(len(pat)):
        if test(i, 0, 0):
            return True
    return False


def verify(result: list[list[int]], solution: list[list[int]]) -> bool:
    if len(result) != len(solution):
        return False
    for lr, ls in zip(result, solution):
        if len(lr) != len(ls):
            return False
        if not all(v1 == v2 for v1, v2 in zip(lr,ls)):
            return False
    return True


def solve(seq: str, patterns: list[str]) -> list[list[int]]:
    l = len(patterns)
    M = [[0 for _ in range(l)] for _ in range(l)]

    for i in range(l):
        for j in range(i, l):
            t = patterns[i]
            u = patterns[j]

            res = canInterweave(seq, t, u)
            M[i][j] = int(res)
            M[j][i] = int(res)

    return M


def load_results(path: str) -> list[list[int]]:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    interwoven_matrix: list[list[int]] = []
    for line in lines:
        vals = [int(v) for v in line.split()]
        interwoven_matrix.append(vals)
    return interwoven_matrix


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq = lines[0]
    patterns = lines[1:]

    M = solve(seq, patterns)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(M, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_itwv_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    patterns = lines[1:]

    M = solve(seq, patterns)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for line in M:
        out = ' '.join(str(int(v)) for v in line)
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
