from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

def parseLine(line: str) -> set[int]:
    res: set[int] = set(int(val) for val in line[1:-1].split(', '))
    return res


"""
https://rosalind.info/problems/seto/

If A and B are sets, then their union A∪B is the set comprising any elements in either A or B; their intersection A∩B is the set of elements in both A and B; and their set difference A−B is the set of elements in A but not in B.

Furthermore, if A is a subset of another set U, then the set complement of A with respect to U is defined as the set Ac=U−A. See the Sample sections below for examples.

    Given: A positive integer n (n≤20,000) and two subsets A and B of {1,2,…,n}.

    Return: Six sets: A∪B, A∩B, A-B, B-A, Ac, and Bc (where set complements are taken with respect to {1,2,…,n}).
"""


def verify(result: list[set[int]], solution: list[set[int]]) -> bool:
    correct = len(result) == len(solution) and all(s1 == s2 for s1, s2 in zip(result, solution))
    return correct


def solve(n: int, set1: set[int], set2: set[int]) -> list[set[int]]:
    completeSet = set(range(1, n+1))
    union = set1.union(set2)
    intersection = set1.intersection(set2)
    diff1 = set1.difference(set2)
    diff2 = set2.difference(set1)
    set1c = completeSet.difference(set1)
    set2c = completeSet.difference(set2)

    return [union, intersection, diff1, diff2, set1c, set2c]


def load_results(path: str) -> list[set[int]]:
    lines = readTextFile(path)
    sets: list[set[int]] = [parseLine(line) for line in lines if len(line) and not line.isspace()]
    return sets


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = int(lines[0])
    set1 = parseLine(lines[1])
    set2 = parseLine(lines[2])

    sets = solve(n, set1, set2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(sets, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_seto_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    set1 = parseLine(lines[1])
    set2 = parseLine(lines[2])

    sets = solve(n, set1, set2)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for set_ in sets:
        print(set_)
        writeTextFile(result_path, str(set_), 'a')

    # correct = solve_and_check(path)
    # print(correct)
