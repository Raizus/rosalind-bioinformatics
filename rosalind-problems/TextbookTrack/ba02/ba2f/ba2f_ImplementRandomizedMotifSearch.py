from BioInfoToolkit.Sequences.SequenceUtils import randomizedMotifSearchMonteCarlo
from BioInfoToolkit.Sequences.StringUtils import getMinimumHammingDistanceToKmer
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba2f/

We will now turn to randomized algorithms that flip coins and roll dice in order to search for motifs. Making random algorithmic decisions may sound like a disastrous idea; just imagine a chess game in which every move would be decided by rolling a die. However, an 18th Century French mathematician and naturalist, Comte de Buffon, first proved that randomized algorithms are useful by randomly dropping needles onto parallel strips of wood and using the results of this experiment to accurately approximate the constant π.

Randomized algorithms may be nonintuitive because they lack the control of traditional algorithms. Some randomized algorithms are Las Vegas algorithms, which deliver solutions that are guaranteed to be exact, despite the fact that they rely on making random decisions. Yet most randomized algorithms are Monte Carlo algorithms. These algorithms are not guaranteed to return exact solutions, but they do quickly find approximate solutions. Because of their speed, they can be run many times, allowing us to choose the best approximation from thousands of runs.

A randomized algorithm for motif finding is given below.

    RANDOMIZEDMOTIFSEARCH(Dna, k, t)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string
            from Dna
        BestMotifs ← Motifs
        while forever
            Profile ← Profile(Motifs)
            Motifs ← Motifs(Profile, Dna)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
            else
                return BestMotifs

Implement RandomizedMotifSearch

    Given: Positive integers k and t, followed by a collection of strings Dna.

    Return: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1000 times. Remember to use pseudocounts!
"""

OutputT = list[str]


def verify(result: OutputT, solution: OutputT, sequences: list[str]) -> bool:
    score_res = sum(sum(getMinimumHammingDistanceToKmer(seq, kmer)
                        for seq in sequences) for kmer in result)
    score_sol = sum(sum(getMinimumHammingDistanceToKmer(seq, kmer)
                        for seq in sequences) for kmer in solution)

    # TODO: change verification:
    # might not be equal since we are doing a monte carlo search

    # alphabet = ['A', 'C', 'G', 'T']
    # profile = SequencesProfile(result, alphabet, 1)
    # score_res_2 = profile.score()

    correct = score_res == score_sol and len(result) == len(solution)
    return correct


def solve(sequences: list[str], k: int) -> OutputT:
    alphabet = ['A', 'C', 'G', 'T']
    n_iter = 1000
    pseudocounts = 1
    best_motifs, _ = randomizedMotifSearchMonteCarlo(
        sequences, k, alphabet, n_iter, pseudocounts)
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
    path = f'{cwd}/rosalind_ba2f_1.txt'

    lines = readTextFile(path)
    k, t = [int(a) for a in lines[0].split()]
    sequences = lines[1:1+t]
    alphabet = ['A', 'C', 'G', 'T']

    best_motifs = solve(sequences, k)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for motif in best_motifs:
        print(motif)
        writeTextFile(result_path, motif, 'a')

    correct = solve_and_check(path)
    print(correct)
