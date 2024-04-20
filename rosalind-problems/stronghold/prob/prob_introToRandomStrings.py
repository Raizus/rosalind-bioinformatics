from math import log
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
An array is a structure containing an ordered collection of objects (numbers, strings, other arrays, etc.). We let A[k] denote the k-th value in array A. You may like to think of an array as simply a matrix having only one row.

A random string is constructed so that the probability of choosing each subsequent symbol is based on a fixed underlying symbol frequency. GC-content offers us natural symbol frequencies for constructing random DNA strings. If the GC-content is x, then we set the symbol frequencies of C and G equal to x2 and the symbol frequencies of A and T equal to 1-x2. For example, if the GC-content is 40%, then as we construct the string, the next symbol is 'G'/'C' with probability 0.2, and the next symbol is 'A'/'T' with probability 0.3.

In practice, many probabilities wind up being very small. In order to work with small probabilities, we may plug them into a function that "blows them up" for the sake of comparison. Specifically, the common logarithm of x
(defined for x>0 and denoted log10(x)) is the exponent to which we must raise 10 to obtain x.

See Figure 1 for a graph of the common logarithm function y=log10(x). In this graph, we can see that the logarithm of x-values between 0 and 1 always winds up mapping to y-values between −∞ and 0: x-values near 0 have logarithms close to −∞, and x-values close to 1 have logarithms close to 0. Thus, we will select the common logarithm as our function to "blow up" small probability values for comparison.

    Given: A DNA string s of length at most 100 bp and an array A containing at most 20 numbers between 0 and 1.

    Return: An array B having the same length as A in which B[k] represents the common logarithm of the probability that a random string constructed with the GC-content found in A[k] will match s exactly.
"""


def verify(result: list[float], solution: list[float]) -> bool:
    correct = all(abs(v1-v2)<=0.001 for v1, v2 in zip(result, solution))
    return correct


def solve(seq: str, Ak: list[float]) -> list[float]:
    probs: list[float] = []
    for ak in Ak:
        prob_dict: dict[str, float] = dict()
        prob_dict["G"] = ak/2
        prob_dict["C"] = ak/2
        prob_dict["A"] = (1-ak)/2
        prob_dict["T"] = (1-ak)/2

        log_likelihood = 0
        for s in seq:
            log_likelihood += log(prob_dict[s], 10)
        probs.append(log_likelihood)
    return probs


def load_results(path: str) -> list[float]:
    lines = readTextFile(path)
    probs = [float(val) for val in lines[0].split(" ")]
    return probs


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    seq = lines[0]
    Ak = [float(val) for val in lines[1].split(" ")]

    result = solve(seq, Ak)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_prob_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    Ak = [float(val) for val in lines[1].split(" ")]

    probs = solve(seq, Ak)

    result_path = result_path_from_input_path(path)

    out = " ".join(f"{p:.4f}" for p in probs)
    print(out)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
