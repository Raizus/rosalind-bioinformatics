import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
    Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

    Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
"""

def prob_dominant(k: int, m: int, n: int) -> float:
    """Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive

    Returns The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
    Args:
        k (int): _description_
        m (int): _description_
        n (int): _description_

    Returns:
        float: _description_
    """
    t = k + m + n

    prob = k/t*((k-1 + m + n)/(t-1)) \
        + m/t * ((k+(m-1)*3/4.0+n/2.0)/(t-1)) \
        + n/t * ((k+m/2.0)/(t-1))
    return prob



def verify(result: float, solution: float) -> bool:
    error = abs(result - solution)
    correct = error <= 0.001
    return correct


def solve(k: int, m: int, n: int) -> float:
    prob = prob_dominant(k, m, n)
    return prob


def load_results(path: str) -> float:
    lines = readTextFile(path)
    prob = float(lines[0])
    return prob


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    (n, k, m) = [int(v) for v in lines[0].split()]

    prob = solve(n, k, m)
    solution = load_results(solution_path)

    correct = verify(prob, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_iprb_1.txt'

    lines = readTextFile(path)
    (n,k,m) = [int(v) for v in lines[0].split()]

    prob = solve(n, k, m)
    print(prob)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(prob), 'w')
