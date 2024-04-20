from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.structures import RNA_CODONS_REVERSE
import os

"""
For positive integers a and n, a modulo n (written amodn in shorthand) is the remainder when a is divided by n. For example, 29mod11=7 because 29=11×2+7.

Modular arithmetic is the study of addition, subtraction, multiplication, and division with respect to the modulo operation. We say that a and b are congruent modulo n if amodn=bmodn; in this case, we use the notation a≡bmodn.

Two useful facts in modular arithmetic are that if a≡bmodn and c≡dmodn, then a+c≡b+dmodn and a×c≡b×dmodn. To check your understanding of these rules, you may wish to verify these relationships for a=29, b=73, c=10, d=32, and n=11.

As you will see in this exercise, some Rosalind problems will ask for a (very large) integer solution modulo a smaller number to avoid the computational pitfalls that arise with storing such large numbers.

    Given: A protein string of length at most 1000 aa.

    Return: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)
"""


def countRNAstringsFromProtein(protein: str, mod: int):
    count = 1
    for aa in protein:
        codons = RNA_CODONS_REVERSE[aa]
        count = (count*len(codons)) % mod
    count = (count*len(RNA_CODONS_REVERSE["_"])) % mod
    return count


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(prot: str) -> int:
    mod = 1000000
    count = countRNAstringsFromProtein(protein, mod)
    return count


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    protein = lines[0]
    result = solve(protein)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_mrna_1.txt'

    mod = 1000000
    lines = readTextFile(path)
    protein = lines[0]

    count = countRNAstringsFromProtein(protein, mod)
    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')
