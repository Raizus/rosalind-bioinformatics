from BioInfoToolkit.Sequences.Reversals import getReversalDistanceWithHistories, performReversals, reversalDistance
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/sort/

A reversal of a permutation can be encoded by the two indices at the endpoints of the interval that it inverts; for example, the reversal that transforms (4,1,2,6,3,5) into (4,1,3,6,2,5) is encoded by [3,5].

A collection of reversals sorts π into γ if the collection contains drev(π,γ) reversals, which when successively applied to π yield γ.

    Given: Two permutations π and γ, each of length 10.

    Return: The reversal distance drev(π,γ), followed by a collection of reversals sorting π into γ. If multiple collections of such reversals exist, you may return any one.
"""


def verify(result: tuple[int, list[tuple[int, int]]], seq1: list[int], seq2: list[int], dist_solution: int) -> bool:
    dist_result = result[0]
    reversals_result = result[1]

    reversed_seq = performReversals(seq1, reversals_result)

    correct = reversed_seq == seq2 and dist_result == dist_solution
    return correct


def solve(perm1: list[int], perm2: list[int]) -> tuple[int, list[tuple[int, int]]]:
    dist, histories = getReversalDistanceWithHistories(perm1, perm2)

    reversals: list[tuple[int, int]] = []

    history = histories[0]
    for rev_start, rev_end in history[1]:
        reversals.append((rev_start, rev_end))

    return dist, reversals


def load_results(path: str) -> tuple[int, list[tuple[int, int]]]:
    lines = readTextFile(path)
    dist = int(lines[0])
    reversals: list[tuple[int, int]] = []

    for line in lines[1:]:
        if not len(line) or line.isspace():
            continue
        i1, i2 = [int(v) for v in line.split()]
        reversals.append((i1-1, i2))

    return dist, reversals


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    perm1 = [int(a) for a in lines[0].split(' ')]
    perm2 = [int(a) for a in lines[1].split(' ')]

    dist, reversals = solve(perm1, perm2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify((dist, reversals), perm1, perm2, solution[0])
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_sort_1.txt'

    lines = readTextFile(path)
    perm1 = [int(a) for a in lines[0].split(' ')]
    perm2 = [int(a) for a in lines[1].split(' ')]

    dist, reversals = solve(perm1, perm2)

    print(dist)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(dist), 'w')

    for rev_start, rev_end in reversals:
        out = f"{rev_start+1} {rev_end}"
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
