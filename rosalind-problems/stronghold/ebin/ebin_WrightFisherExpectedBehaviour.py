import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ebin/

In “The Wright-Fisher Model of Genetic Drift”, we generalized the concept of a binomial random variable Bin(n,p) as a "weighted coin flip." It is only natural to calculate the expected value of such a random variable.

For example, in the case of unweighted coin flips (i.e., p=1/2), our intuition would indicate that E(Bin(n,1/2)) is n/2; what should be the expected value of a binomial random variable?

    Given: A positive integer n (n≤1000000) followed by an array P of length m (m≤20) containing numbers between 0 and 1. Each element of P can be seen as representing a probability corresponding to an allele frequency.

    Return: An array B of length m for which B[k] is the expected value of Bin(n,P[k]); in terms of Wright-Fisher, it represents the expected allele frequency of the next generation.
"""

OutputT = list[float]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = len(result) == len(solution) and all(abs(v1-v2) <= 0.001 for v1, v2 in zip(result, solution))
    return correct


def solve(n: int, P: list[float]) -> OutputT:
    mean = [n*ak for ak in P]
    return mean


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    Bk = [float(v) for v in lines[0].split()]
    return Bk


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = int(lines[0])
    P = [float(v) for v in lines[1].split()]

    result = solve(n, P)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ebin_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    P = [float(v) for v in lines[1].split()]

    mean = solve(n, P)

    out = ' '.join(f"{v:.4f}" for v in mean)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

