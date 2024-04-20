from BioInfoToolkit.IO import readTextFile, writeTextFile
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/lexv/

Say that we have strings s=s1s2⋯sm and t=t1t2⋯tn with m<n. Consider the substring t′=t[1:m]. We have two cases:

    1. If s=t', then we set s<Lext because s is shorter than t (e.g., APPLE<APPLET).
    2. Otherwise, s≠t'. We define s<Lext if s<Lext' and define s>Lext if s>Lext' (e.g., APPLET<LexARTS because APPL<LexARTS).

    Given: A permutation of at most 12 symbols defining an ordered alphabet A and a positive integer n (n≤4).

    Return: All strings of length at most n formed from A, ordered lexicographically. (Note: As in “Enumerating k-mers Lexicographically”, alphabet order is based on the order in which the symbols are given.)
"""


def generateOrderedString(alphabet: list[str], max_length: int, prev: str = ""):
    for s in alphabet:
        res = prev + s
        yield res
        if len(res) < max_length:
            for res2 in generateOrderedString(alphabet, max_length, res):
                yield res2


def verify(result: list[str], solution: list[str]) -> bool:
    correct = len(result) == len(solution) and all(w1 == w2 for w1, w2 in zip(result, solution))
    return correct


def solve(alphabet: list[str], n: int) -> list[str]:
    words: list[str] = []
    for string in generateOrderedString(alphabet, n, ""):
        words.append(string)
    return words


def load_results(path: str) -> list[str]:
    words = readTextFile(path)
    words = [line for line in words if len(line) and not line.isspace()]
    return words


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
    path = f'{cwd}/rosalind_lexv_1.txt'

    lines = readTextFile(path)
    alphabet = lines[0].split()
    n = int(lines[1])

    words = solve(alphabet, n)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for word in words:
        print(word)
        writeTextFile(result_path, word, 'a')

    # correct = solve_and_check(path)
    # print(correct)
