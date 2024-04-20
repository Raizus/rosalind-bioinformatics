
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from scipy.stats import binom
import numpy as np

"""
https://rosalind.info/problems/wfmd/

Consider flipping a weighted coin that gives "heads" with some fixed probability p (i.e., p is not necessarily equal to 1/2).

We generalize the notion of binomial random variable from “Independent Segregation of Chromosomes” to quantify the sum of the weighted coin flips. Such a random variable X takes a value of k if a sequence of n independent "weighted coin flips" yields k "heads" and n-k "tails." We write that X∈Bin(n,p).

To quantify the Wright-Fisher Model of genetic drift, consider a population of N diploid individuals, whose 2N chromosomes possess m copies of the dominant allele. As in “Counting Disease Carriers”, set p=m2N. Next, recall that the next generation must contain exactly N individuals. These individuals' 2N alleles are selected independently: a dominant allele is chosen with probability p, and a recessive allele is chosen with probability 1-p.

    Given: Positive integers N (N≤7), m (m≤2N), g (g≤6) and k (k≤2N).

    Return: The probability that in a population of N diploid individuals initially possessing m copies of a dominant allele, we will observe after g generations at least k copies of a recessive allele. Assume the Wright-Fisher model.
"""


def pmf_scipy(n, p):
    pmfA = [float(binom.pmf(k, n, p)) for k in range(0, n+1)]
    return pmfA


def cdf_scipy(n, p):
    cdfA = [float(binom.cdf(k, n, p)) for k in range(0, n+1)]
    return cdfA


def getTransitionMatrix(n: int):
    P = np.zeros((n+1, n+1))
    for i in range(n+1):
        for j in range(n+1):
            P[i][j] = float(binom.pmf(j, n, i/n))
    return P


def wrightFisherNextGen(p: list[float], n: int):
    assert len(p) == n+1
    return [sum([float(binom.pmf(r, n, i/n))*p[i] for i in range(n+1)])
            for r in range(1+n)]


def wrightFisherModel(n: int, m: int, g: int) -> list[float]:
    """Cumputes the probability mass function pmf(x), where P(X=k), k in [0,m], is the probability that there are k alleles in the population after g generations

    Args:
        n (int): Number of chromosomes in the population
        m (int): Number of alleles in the population
        g (int): Number of gerenations

    Returns:
        list[float]: _description_
    """
    ms: list[float] = [1.0 if i == m else 0 for i in range(n+1)]
    for _ in range(g):
        ms = wrightFisherNextGen(ms, n)
    return ms


def verify(result: float, solution: float) -> bool:
    correct = abs(result-solution) <= 0.001
    return correct


def solve(N: int, m: int, g: int, k: int) -> float:
    totalAlleles = 2*N
    m_recessive = totalAlleles - m

    # p_recessive_allele = m_recessive/totalAlleles
    # A = np.zeros((totalAlleles + 1, totalAlleles + 1))
    # A[m_recessive][m_recessive] = 1

    # P = getTransitionMatrix(totalAlleles)
    # res = A @ np.linalg.matrix_power(P, g)
    # prob = sum(res[m_recessive][k:])

    ms = wrightFisherModel(2*N, m_recessive, g)
    prob = sum(ms[k:])

    return prob


def load_results(path: str) -> float:
    lines = readTextFile(path)
    prob = float(lines[0])
    return prob


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n, m, g, k = [int(v) for v in lines[0].split()]

    prob = solve(n, m, g, k)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(prob, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_wfmd_1.txt'

    lines = readTextFile(path)
    n, m, g, k = [int(v) for v in lines[0].split()]

    prob = solve(n, m, g, k)

    print(prob)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(prob), 'w')

    correct = solve_and_check(path)
    print(correct)

