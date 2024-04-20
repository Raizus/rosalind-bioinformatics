from BioInfoToolkit.Sequences.SequenceUtils import gibbsSamplerMonteCarlo
from BioInfoToolkit.Sequences.StringUtils import getMinimumHammingDistanceToKmer
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba2g/

We have previously defined the notion of a Profile-most probable k-mer in a string. We now define a Profile-randomly generated k-mer in a string Text. For each k-mer Pattern in Text, compute the probability Pr(Pattern | Profile), resulting in n = |Text| - k + 1 probabilities (p1, …, pn). These probabilities do not necessarily sum to 1, but we can still form the random number generator Random(p1, …, pn) based on them. GIBBSSAMPLER uses this random number generator to select a Profile-randomly generated k-mer at each step: if the die rolls the number i, then we define the Profile-randomly generated k-mer as the i-th k-mer in Text.

    GIBBSSAMPLER(Dna, k, t, N)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string
            from Dna
        BestMotifs ← Motifs
        for j ← 1 to N
            i ← Random(t)
            Profile ← profile matrix constructed from all strings in Motifs
                       except for Motifi
            Motifi ← Profile-randomly generated k-mer in the i-th sequence
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs

Implement GibbsSampler

    Given: Integers k, t, and N, followed by a collection of strings Dna.

    Return: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!

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

    correct = score_res <= score_sol and len(result) == len(solution)
    return correct


def solve(sequences: list[str], k: int, N: int) -> OutputT:
    alphabet = ['A', 'C', 'G', 'T']
    pseudocounts = 1
    n_starts = 20

    best_motifs, _ = gibbsSamplerMonteCarlo(
        sequences, k, alphabet, N, n_starts, pseudocounts)
    return best_motifs


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    motifs = [line for line in lines if len(line) and not line.isspace()]
    return motifs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    k, t, N = [int(a) for a in lines[0].split()]
    sequences = lines[1:1+t]

    result = solve(sequences, k, N)

    solution = load_results(solution_path)

    correct = verify(result, solution, sequences)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba2g_1.txt'

    lines = readTextFile(path)
    k, t, N = [int(a) for a in lines[0].split()]
    sequences = lines[1:1+t]
    alphabet = ['A', 'C', 'G', 'T']

    best_motifs = solve(sequences, k, N)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for motif in best_motifs:
        print(motif)
        writeTextFile(result_path, motif, 'a')

    correct = solve_and_check(path)
    print(correct)
