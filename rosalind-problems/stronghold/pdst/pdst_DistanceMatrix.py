from BioInfoToolkit.Sequences.SequenceUtils import distanceMatrix
import os
from BioInfoToolkit.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/pdst/

For two strings s1 and s2 of equal length, the p-distance between them, denoted dp(s1,s2), is the proportion of corresponding symbols that differ between s1 and s2.

For a general distance function d on n taxa s1,s2,…,sn (taxa are often represented by genetic strings), we may encode the distances between pairs of taxa via a distance matrix D in which Di,j=d(si,sj).

    Given: A collection of n (n≤10) DNA strings s1,…,sn of equal length (at most 1 kbp). Strings are given in FASTA format.

    Return: The matrix D corresponding to the p-distance dp on the given strings. As always, note that your answer is allowed an absolute error of 0.001.
"""


def verify(result: list[list[float]], solution: list[list[float]]) -> bool:
    if len(result) != len(solution):
        return False
    for line1, line2 in zip(result, solution):
        if len(line1) != len(line2) or not all(abs(v1-v2) <= 0.001 for v1, v2 in zip(line1, line2)):
            return False
    return True


def solve(seqs: list[str]) -> list[list[float]]:
    dist_mat = distanceMatrix(seqs)
    dist_mat = dist_mat.tolist()
    return dist_mat


def load_results(path: str) -> list[list[float]]:
    lines = readTextFile(path)
    dist_mat: list[list[float]] = [[float(v) for v in line.split()] for line in lines]
    return dist_mat


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    dist_mat = solve(seqs)

    solution = load_results(solution_path)
    correct = verify(dist_mat, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_pdst_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = list(fasta_dict.values())
    dist_mat = solve(seqs)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for row in dist_mat:
        out = " ".join(f"{v:.4f}" for v in row)
        print(out)
        writeTextFile(result_path, out, 'a')

    # correct = solve_and_check(path)
    # print(correct)

