from BioInfoToolkit.Sequences.BioSequence import BioSequence
from BioInfoToolkit.IO import read_FASTA

from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
A DNA string is a reverse palindrome if it is equal to its reverse complement. For instance, GCATGC is a reverse palindrome because its reverse complement is GCATGC.

    Given: A DNA string of length at most 1 kbp in FASTA format.

    Return: The position and length of every reverse palindrome in the string having length between 4 and 12. You may return these pairs in any order.
"""


def verify(result: list[tuple[int, int]], solution: list[tuple[int, int]]) -> bool:
    for pair in result:
        if pair not in solution:
            return False

    return True


def solve(seq: str, la: int, lb: int) -> list[tuple[int, int]]:
    pos_len_pairs: list[tuple[int, int]] = []
    bio_seq = BioSequence(seq)

    for l in range(la, lb+1):
        res = bio_seq.find_reverse_palindrome(l)
        for i in res:
            pos_len_pairs.append((i+1, l))
    return pos_len_pairs


def load_results(path: str) -> list[tuple[int, int]]:
    lines = readTextFile(path)
    pos_len_pairs: list[tuple[int, int]] = []
    for line in lines:
        if not len(line):
            continue
        p, l = line.split()
        pos_len_pairs.append((int(p), int(l)))
    # pos_len_pairs = [(int(pl[0]), int(pl[1])) for line in lines for pl in line.split() if len(line)]
    return pos_len_pairs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    la, lb = 4, 12
    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]
    pos_len_pairs = solve(seq, la, lb)

    solution = load_results(solution_path)

    correct = verify(pos_len_pairs, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_revp_1.txt'

    la, lb = 4, 12
    
    fasta_dict = read_FASTA(path)
    seq = list(fasta_dict.values())[0]
    pos_len_pairs = solve(seq, la, lb)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for (p, l) in pos_len_pairs:
        out = f"{p} {l}"
        print(out)
        writeTextFile(result_path, out, 'a')
