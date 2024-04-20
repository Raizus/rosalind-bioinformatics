from BioInfoToolkit.Sequences.Reversals import reversalDistance
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
A reversal of a permutation creates a new permutation by inverting some interval of the permutation; (5,2,3,1,4), (5,3,4,1,2), and (4,1,2,3,5) are all reversals of (5,3,2,1,4). The reversal distance between two permutations π and sigma, written drev(π,sigma), is the minimum number of reversals required to transform π into sigma (this assumes that π and sigma have the same length).

    Given: A collection of at most 5 pairs of permutations, all of which have length 10.

    Return: The reversal distance between each permutation pair.
"""


def verify(result: list[int], solution: list[int]) -> bool:
    correct = result == solution
    return correct


def solve(perm_pairs: list[tuple[list[int], list[int]]]) -> list[int]:
    distances: list[int] = []
    for perm1, perm2 in perm_pairs:
        dist = reversalDistance(perm1, perm2)
        distances.append(dist)
    return distances


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    distances = [int(v) for v in lines[0].split()]
    return distances


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = len(lines)
    perm_pairs: list[tuple[list[int], list[int]]] = []

    for i in range(0, n, 3):
        seq1 = [int(val) for val in lines[i].split(' ')]
        seq2 = [int(val) for val in lines[i+1].split(' ')]
        pair = (seq1, seq2)
        perm_pairs.append(pair)

    distances = solve(perm_pairs)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(distances, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_rear_0.txt'

    # TODO Fix This

    lines = readTextFile(path)
    n = len(lines)
    perm_pairs: list[tuple[list[int], list[int]]] = []

    for i in range(0, n, 3):
        seq1 = [int(val) for val in lines[i].split(' ')]
        seq2 = [int(val) for val in lines[i+1].split(' ')]
        pair = (seq1, seq2)
        perm_pairs.append(pair)

    distances = solve(perm_pairs)

    out = ' '.join(str(d) for d in distances)

    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)