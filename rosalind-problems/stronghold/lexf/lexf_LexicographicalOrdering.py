
from itertools import product
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
In a weighted alphabet, every symbol is assigned a positive real number called a weight. A string formed from a weighted alphabet is called a weighted string, and its weight is equal to the sum of the weights of its symbols.

The standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic mass of the corresponding amino acid.

    Given: A protein string P of length at most 1000 aa.

    Return: The total weight of P. Consult the monoisotopic mass table.
"""


def verify(result: list[str], solution: list[str]) -> bool:
    correct = all(w1 == w2 for w1, w2 in zip(result, solution))
    return correct


def solve(alphabet: list[str], n: int) -> list[str]:
    words: list[str] = []
    for p in product(alphabet, repeat=n):
        w = "".join(p)
        words.append(w)
    return words


def load_results(path: str) -> list[str]:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    return lines


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    alphabet = lines[0].split()
    n = int(lines[1])

    words = solve(alphabet, n)
    solution = load_results(solution_path)

    correct = verify(words, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_lexf_1.txt'

    lines = readTextFile(path)
    alphabet = lines[0].split()
    n = int(lines[1])

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for p in product(alphabet, repeat=n):
        out = "".join(p)
        print(out)
        writeTextFile(result_path, out, 'a')

    # correct = solve_and_check(path)
    # print(correct)