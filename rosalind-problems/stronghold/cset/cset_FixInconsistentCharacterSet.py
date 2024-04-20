from itertools import combinations
from BioInfoToolkit.Phylogeny.Phylo import conflicting_splits, is_character_table_consistent

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/cset/

A submatrix of a matrix M is a matrix formed by selecting rows and columns from M and taking only those entries found at the intersections of the selected rows and columns. We may also think of a submatrix as formed by deleting the remaining rows and columns from M.

Given: An inconsistent character table C on at most 100 taxa.

Return: A submatrix of C' representing a consistent character table on the same taxa and formed by deleting a single row of C. (If multiple solutions exist, you may return any one.)
"""


def remove_inconsistent_line(character_table: list[str]):
    character_table = character_table.copy()
    ntaxa = len(character_table[0])
    nchars = len(character_table)

    for i1, i2 in combinations(range(nchars), 2):
        ch1 = character_table[i1]
        ch2 = character_table[i2]

        S1 = set(i for i, c in enumerate(ch1) if c == '1')
        S1c = set(range(ntaxa)).difference(S1)
        S2 = set(i for i, c in enumerate(ch2) if c == '1')
        S2c = set(range(ntaxa)).difference(S2)
        inconsistent = conflicting_splits(S1, S1c, S2, S2c)

        if inconsistent:
            character_table = character_table[:i2] + character_table[i2+1:]
            break

    return character_table


OutputT = list[str]


def verify(result: OutputT, solution: OutputT, input: list[str]) -> bool:
    correct = is_character_table_consistent(result) and \
        len(result) == len(input) - 1 and \
        len(set(result).intersection(set(input))) == len(input) - 1
    return correct


def solve(character_table: list[str]) -> OutputT:
    consistent_character_table = remove_inconsistent_line(character_table)
    return consistent_character_table


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    result = [line for line in lines if len(line) and not line.isspace()]
    return result


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    character_table = [line for line in lines if len(
        line) and not line.isspace()]

    result = solve(character_table)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, character_table)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_cset_1.txt'

    lines = readTextFile(path)
    character_table = [line for line in lines if len(
        line) and not line.isspace()]

    # consistency_matrix(character_table)
    consistent_character_table = remove_inconsistent_line(character_table)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for c in consistent_character_table:
        print(c)
        writeTextFile(result_path, c, 'a')

    correct = solve_and_check(path)
    print(correct)
