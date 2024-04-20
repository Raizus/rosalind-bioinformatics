import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/rstr/

Our aim in this problem is to determine the probability with which a given motif (a known promoter, say) occurs in a randomly constructed genome. Unfortunately, finding this probability is tricky; instead of forming a long genome, we will form a large collection of smaller random strings having the same length as the motif; these smaller strings represent the genome's substrings, which we can then test against our motif.

Given a probabilistic event A, the complement of A is the collection Ac of outcomes not belonging to A. Because Ac takes place precisely when A does not, we may also call Ac "not A."

For a simple example, if A is the event that a rolled die is 2 or 4, then Pr(A)=13. Ac is the event that the die is 1, 3, 5, or 6, and Pr(Ac)=23. In general, for any event we will have the identity that Pr(A)+Pr(Ac)=1.

    Given: A positive integer N≤100000, a number x between 0 and 1, and a DNA string s of length at most 10 bp.

    Return: The probability that if N random DNA strings having the same length as s are constructed with GC-content x (see “Introduction to Random Strings”), then at least one of the strings equals s. We allow for the same random string to be created more than once.
"""


def verify(result: float, solution: float) -> bool:
    correct = abs(result-solution) <= 0.001
    return correct


def solve(N: int, x: float, seq: str) -> float:
    l = len(seq)
    # First, consider the event E, that a random sequence matches a given motif with length L. With GC content given by x, the vector of probability distribution of A,C,G,T is given by ((1-x)/2, x/2, x/2, (1-x)/2)

    # prob of a random string of length l having gc-content x
    gc_count = seq.count('G') + seq.count('C')
    p_single = (x/2)**gc_count * ((1-x)/2)**(l-gc_count)

    # Compute probability that none of the string equals s P(X = 0)
    p_none_equals_s = (1-p_single)**N
    p_at_least_one_equals_s = 1 - p_none_equals_s
    return p_at_least_one_equals_s


def load_results(path: str) -> float:
    lines = readTextFile(path)
    prob = float(lines[0])
    return prob


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    vals = lines[0].split()
    N = int(vals[0])
    x = float(vals[1])
    seq = lines[1]

    prob = solve(N, x, seq)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(prob, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_rstr_1.txt'

    lines = readTextFile(path)
    vals = lines[0].split()
    N = int(vals[0])
    x = float(vals[1])
    seq = lines[1]

    prob = solve(N,x, seq)

    print(prob)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(prob), 'w')

    # correct = solve_and_check(path)
    # print(correct)
