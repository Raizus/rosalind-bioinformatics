from scipy.stats import binom
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation F_n = F_(n-1)+F_(n-2) and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months.

    Given: Positive integers n≤100 and m≤20.

    Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.
"""


def verify(result: float, solution: float) -> bool:
    error = abs(result - solution)
    correct = error <= 0.001
    return correct



def solve(k: int, N: int, n0: int) -> float:
    """Computes the probability that at least N Aa Bb organisms will belong to the k-th generation of Tom's family tree (don't count the Aa Bb mates at each level). Assume that Mendel's second law holds for the factors.

    Args:
        k (int): number of generations
        N (int): number of Aa Bb organisms
        n0 (int): number of offspring per individual

    Returns:
        float: The probability that at least N Aa Bb organisms will belong to the k-th generation of Tom's family tree (don't count the Aa Bb mates at each level).
    """
    # we want P(Xk >= N), where Xk is the number of organisms having genotype Aa Bb, in generation k
    # start at k = 0, with C0 = 1 organism

    # Each organism always mates with an organism having genotype Aa Bb.
    # Offspring probabilities are as follows (for 1 single offspring)

    # If both parents are Aa Bb

    #   AA BB - 1/4 * 1/4 = 1/16
    #   AA Bb - 1/4 * 1/2 = 1/8 = 2/16
    #   AA bb - 1/4 * 1/4 = 1/16

    #   Aa BB - 1/2 * 1/4 = 1/8 = 2/16
    #   Aa Bb - 1/2 * 1/2 = 1/4 = 4/16 <-- only interested in this probability
    #   Aa bb - 1/2 * 1/4 = 1/8 = 2/16

    #   aa BB - 1/4 * 1/4 = 1/16
    #   aa Bb - 1/4 * 1/2 = 1/8 = 2/16
    #   aa bb - 1/4 * 1/4 = 1/16

    p_Aa_Bb = 1/2 * 1/2  # this probability is the same in each generation
    # At the k'th gen we have 2*k offspring
    # P(Xk >= N) = 1 - P(Xk < N) Xk is the random variable
    nTrials = 2**k
    pf = 1 - float(binom.cdf(N-1, nTrials, p_Aa_Bb))
    return pf


def load_results(path: str) -> float:
    lines = readTextFile(path)
    pf = float(lines[0])
    return pf


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    (k, N) = [int(v) for v in lines[0].split()]
    n0 = 2

    result = solve(k, N, n0)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_lia_1.txt'

    lines = readTextFile(path)
    (k, N) = [int(v) for v in lines[0].split()]
    n0 = 2

    pf = solve(k, N, n0)
    print(pf)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(pf), 'w')
