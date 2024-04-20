from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, writeTextFile
import os
from math import floor

"""
A subsequence of a permutation is a collection of elements of the permutation in the order that they appear. For example, (5, 3, 4) is a subsequence of (5, 1, 3, 4, 2).

A subsequence is increasing if the elements of the subsequence increase, and decreasing if the elements decrease. For example, given the permutation (8, 2, 1, 6, 5, 7, 4, 3, 9), an increasing subsequence is (2, 6, 7, 9), and a decreasing subsequence is (8, 6, 5, 4, 3). You may verify that these two subsequences are as long as possible.

    Given: A positive integer n≤10000 followed by a permutation π of length n.

    Return: A longest increasing subsequence of π, followed by a longest decreasing subsequence of π.
"""


def longestDecreasingSubsequence(X: list[int]) -> list[int]:
    if not X:
        return X

    n = len(X)
    # stores the index k of the smallest value X[k] such that there is an increasing
    # subsequence of length l ending at X[k] in the range k ≤ i.
    # i1 < i2 < ... < il = k such that X[i1] > X[i2] > ... > X[k]
    M: list[int] = [-1] * (n + 1)
    # P[k] — stores the index of the predecessor of X[k] in the longest decreasing subsequence ending at X[k]
    P: list[int] = [-1] * n

    L = 0
    M[0] = -1

    for i, Xi in enumerate(X):
        # Binary search: we want the smallest l <= L
        # such that X[M[l]] > X[i] (default j = 0),
        lower = 1
        upper = L + 1

        while lower < upper:
            mid = lower + floor((upper-lower)/2)  # lower <= mid < upper
            if X[M[mid]] < Xi:
                upper = mid
            else:  # X[M[mid]] >= X[i]
                lower = mid + 1

        newL = lower
        P[i] = M[newL - 1]
        M[newL] = i

        # found longest subsequence yet
        if newL > L:
            L = newL

    S: list[int] = []
    k = M[L]
    for i in range(L-1, -1, -1):
        S.append(X[k])
        k = P[k]
    S.reverse()

    return S


def longestIncreasingSubsequence(X: list[int]) -> list[int]:
    if not X:
        return X

    n = len(X)
    # stores the index k of the smallest value X[k] such that there is an increasing
    # subsequence of length l ending at X[k] in the range k ≤ i.
    # i1 < i2 < ... < il = k such that X[i1] < X[i2] < ... < X[k]
    M: list[int] = [-1] * (n + 1)
    # P[k] — stores the index of the predecessor of X[k] in the longest increasing subsequence ending at X[k]
    P: list[int] = [-1] * n

    L = 0
    M[0] = -1

    for i, Xi in enumerate(X):
        # Binary search: we want the smallest l <= L
        # such that X[M[l]] < X[i] (default j = 0),
        lower = 1
        upper = L + 1

        while lower < upper:
            mid = lower + floor((upper-lower)/2)  # lower <= mid < upper
            if X[M[mid]] >= Xi:
                upper = mid
            else:  # X[M[mid]] < X[i]
                lower = mid + 1

        newL = lower
        P[i] = M[newL - 1]
        M[newL] = i

        # found longest subsequence yet
        if newL > L:
            L = newL

    S: list[int] = []
    k = M[L]
    for i in range(L-1, -1, -1):
        S.append(X[k])
        k = P[k]
    S.reverse()

    return S


def is_subsequence(seq: list[int], subseq: list[int]):
    i1 = 0
    v1 = subseq[i1]
    l = len(subseq)
    for v2 in seq:
        if v2 != v1:
            continue
        i1 += 1
        if i1 == l:
            return True
        v1 = subseq[i1]
    return False



def is_longest_decreasing_subsequence(seq: list[int], subseq: list[int]) -> bool:
    is_decreasing = all(v1 > v2 for v1, v2 in zip(subseq, subseq[1:]))
    if not is_decreasing:
        return False

    if not is_subsequence(seq, subseq):
        return False

    lds = longestDecreasingSubsequence(seq)
    if len(lds) != len(subseq):
        return False

    return True


def is_longest_increasing_subsequence(seq: list[int], subseq: list[int]) -> bool:
    is_increasing = all(v1 < v2 for v1, v2 in zip(subseq, subseq[1:]))
    if not is_increasing:
        return False

    if not is_subsequence(seq, subseq):
        return False

    lis = longestIncreasingSubsequence(seq)
    if len(lis) != len(subseq):
        return False

    return True


def verify(seq: list[int], lis: list[int], lds: list[int]) -> bool:
    a1 = is_longest_increasing_subsequence(seq, lis)
    a2 = is_longest_decreasing_subsequence(seq, lds)
    return a1 and a2


def solve(seq: list[int]) -> tuple[list[int], list[int]]:
    lis = longestIncreasingSubsequence(seq)
    lds = longestDecreasingSubsequence(seq)

    return (lis, lds)


def load_results(path: str) -> tuple[list[int], list[int]]:
    lines = readTextFile(path)
    lis = [int(v) for v in lines[0]]
    lds = [int(v) for v in lines[1]]
    return (lis, lds)


def solve_and_check(input_path: str) -> bool:
    # solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    # n = int(lines[0])
    seq = [int(v) for v in lines[1].split()]

    lis, lds = solve(seq)
    correct = verify(seq, lis, lds)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_lgis_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    seq = [int(v) for v in lines[1].split()]

    lis = longestIncreasingSubsequence(seq)
    lds = longestDecreasingSubsequence(seq)

    lis_out = ' '.join(str(val) for val in lis)
    lds_out = ' '.join(str(val) for val in lds)

    print(lis_out)
    print(lds_out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, lis_out, 'w')
    writeTextFile(result_path, lds_out, 'a')

    correct = solve_and_check(path)
    print(correct)

