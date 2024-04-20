import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/sexl/

The conditional probability of an event A given another event B, written Pr(A|B), is equal to Pr(A and B) divided by Pr(B).

Note that if A and B are independent, then Pr(A and B) must be equal to Pr(A) x Pr(B), which results in Pr(A|B)=Pr(A). This equation offers an intuitive view of independence: the probability of A, given the occurrence of event B, is simply the probability of A (which does not depend on B).

In the context of sex-linked traits, genetic equilibrium requires that the alleles for a gene k are uniformly distributed over the males and females of a population. In other words, the distribution of alleles is independent of sex.

    Given: An array A of length n for which A[k] represents the proportion of males in a population exhibiting the k-th of n total recessive X-linked genes. Assume that the population is in genetic equilibrium for all n genes.

    Return: An array B of length n in which B[k] equals the probability that a randomly selected female will be a carrier for the k-th gene.
"""


def verify(result: list[float], solution: list[float]) -> bool:
    correct = len(result) == len(solution) and all(
        abs(v1-v2) <= 0.001 for v1, v2 in zip(result, solution))
    return correct


def solve(Ak: list[float]) -> list[float]:
    Bk = [2*ak*(1-ak) for ak in Ak]
    return Bk


def load_results(path: str) -> list[float]:
    lines = readTextFile(path)
    Bk = [float(v) for v in lines[0].split()]
    return Bk


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    Ak = [float(v) for v in lines[0].split()]

    probs = solve(Ak)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(probs, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_sexl_1.txt'

    lines = readTextFile(path)
    Ak = [float(v) for v in lines[0].split()]

    probs = solve(Ak)

    out = ' '.join(f"{val:.4f}" for val in probs)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(out), 'w')

    correct = solve_and_check(path)
    print(correct)
