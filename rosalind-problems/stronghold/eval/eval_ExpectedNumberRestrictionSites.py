import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/eval/

Say that you place a number of bets on your favorite sports teams. If their chances of winning are 0.3, 0.8, and 0.6, then you should expect on average to win 0.3 + 0.8 + 0.6 = 1.7 of your bets (of course, you can never win exactly 1.7!)

More generally, if we have a collection of events A1,A2,…,An, then the expected number of events occurring is Pr(A1)+Pr(A2)+⋯+Pr(An) (consult the note following the problem for a precise explanation of this fact). In this problem, we extend the idea of finding an expected number of events to finding the expected number of times that a given string occurs as a substring of a random string.

    Given: A positive integer n (n≤1,000,000), a DNA string s of even length at most 10, and an array A of length at most 20, containing numbers between 0 and 1.

    Return: An array B having the same length as A in which B[i] represents the expected number of times that s will appear as a substring of a random DNA string t of length n, where t is formed with GC-content A[i] (see “Introduction to Random Strings”).
"""


def verify(result: list[float], solution: list[float]) -> bool:
    correct = len(result) == len(solution) and all(abs(v1-v2) <= 0.001 for v1, v2 in zip(result, solution))
    return correct


def solve(n: int, seq: str, Ak: list[float]) -> list[float]:
    l = len(seq)
    gc_count = seq.count('G') + seq.count('C')
    probs: list[float] = []
    for ak in Ak:
        Pe = (ak/2)**gc_count * ((1-ak)/2)**(l-gc_count)
        prob = (n-1)*Pe
        probs.append(prob)
    return probs


def load_results(path: str) -> list[float]:
    lines = readTextFile(path)
    probs = [float(v) for v in lines[0].split()]
    return probs


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = int(lines[0])
    seq = lines[1]
    Ak = [float(v) for v in lines[2].split()]

    probs = solve(n, seq, Ak)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(probs, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_eval_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    seq = lines[1]
    Ak = [float(v) for v in lines[2].split()]

    probs = solve(n, seq, Ak)

    out = ' '.join(f"{val:.4f}" for val in probs)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(out), 'w')

    # correct = solve_and_check(path)
    # print(correct)
