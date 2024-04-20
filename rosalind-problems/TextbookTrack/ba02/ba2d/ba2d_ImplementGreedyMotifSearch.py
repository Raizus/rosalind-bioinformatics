from BioInfoToolkit.Sequences.SequenceUtils import greedyMotifSearch
from BioInfoToolkit.Sequences.StringUtils import getMinimumHammingDistanceToKmer
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba2d/

GREEDYMOTIFSEARCH(Dna, k, t)
        BestMotifs ← motif matrix formed by first k-mers in each string
                      from Dna
        for each k-mer Motif in the first string from Dna
            Motif1 ← Motif
            for i = 2 to t
                form Profile from motifs Motif1, …, Motifi - 1
                Motifi ← Profile-most probable k-mer in the i-th string
                          in Dna
            Motifs ← (Motif1, …, Motift)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs

Implement GreedyMotifSearch

    Given: Integers k and t, followed by a collection of strings Dna.

    Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t). If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT, sequences: list[str]) -> bool:
    score_res = sum(sum(getMinimumHammingDistanceToKmer(seq, kmer)
        for seq in sequences) for kmer in result)
    score_sol = sum(sum(getMinimumHammingDistanceToKmer(seq, kmer)
        for seq in sequences) for kmer in solution)

    correct = score_res == score_sol and len(result) == len(solution)
    return correct


def solve(sequences: list[str], k: int) -> OutputT:
    alphabet = ['A', 'C', 'G', 'T']
    best_motifs, _ = greedyMotifSearch(sequences, k, alphabet)
    return best_motifs


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    motifs = [line for line in lines if len(line) and not line.isspace()]
    return motifs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    k, t = [int(a) for a in lines[0].split()]
    sequences = lines[1:1+t]

    result = solve(sequences, k)

    solution = load_results(solution_path)

    correct = verify(result, solution, sequences)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba2d_1.txt'

    lines = readTextFile(path)
    k, t = [int(a) for a in lines[0].split()]
    sequences = lines[1:1+t]

    best_motifs = solve(sequences, k)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for motif in best_motifs:
        print(motif)
        writeTextFile(result_path, motif, 'a')

    correct = solve_and_check(path)
    print(correct)

