import math
from scipy.stats import binom
from decimal import Decimal

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/indc/

Consider a collection of coin flips. One of the most natural questions we can ask is if we flip a coin 92 times, what is the probability of obtaining 51 "heads", vs. 27 "heads", vs. 92 "heads"?

Each coin flip can be modeled by a uniform random variable in which each of the two outcomes ("heads" and "tails") has probability equal to 1/2. We may assume that these random variables are independent (see “Independent Alleles”); in layman's terms, the outcomes of the two coin flips do not influence each other.

A binomial random variable X takes a value of k if n consecutive "coin flips" result in k total "heads" and n-k total "tails." We write that X∈Bin(n,1/2).

    Given: A positive integer n≤50.

    Return: An array A of length 2n in which A[k] represents the common logarithm of the probability that two diploid siblings share at least k of their 2n chromosomes (we do not consider recombination for now).
"""


def pmf_decimal(k, n, p):
    return Decimal(math.comb(n, k))*Decimal(p**k)*Decimal((1-p)**(n-k))


def get_cdf_decimal(n, p):
    pmfA = [pmf_decimal(k, n, p) for k in range(0, n+1)]

    prev = Decimal(0)
    cdf1: list[Decimal] = []
    for i, k in enumerate(range(0, n+1)):
        res = pmfA[i] + prev
        cdf1.append(res)
        prev = res

    cdf2: list[Decimal] = [(val).log10() for val in cdf1[n-1::-1]]
    return cdf2


def get_cdf_scipy(n, p):
    cdfA = [float(binom.cdf(k, n, p)) for k in range(0, n+1)]
    cdfA = [math.log10(val) for val in cdfA]
    return cdfA[n-1::-1]


def cdf(k, n, p):
    return binom.cdf(k, n, p)



def verify(result: list[float], solution: list[float]) -> bool:
    correct = len(result) == len(solution) and all(abs(v1-v2)<=0.001 for v1, v2 in zip(result, solution))
    return correct


def solve(n: int, p: float) -> list[float]:
    cdf_decimal = get_cdf_decimal(2*n, p)
    cdf_scipy = get_cdf_scipy(2*n, p)

    probs = [float(v) for v in cdf_decimal]
    return probs


def load_results(path: str) -> list[float]:
    lines = readTextFile(path)
    probs = [float(v) for v in lines[0].split()]

    return probs


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = int(lines[0])
    p = 0.5

    probs = solve(n, p)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(probs, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_indc_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    p = 0.5

    probs = solve(n, p)
    out = " ".join(f"{v:.4f}" for v in probs)

    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

    