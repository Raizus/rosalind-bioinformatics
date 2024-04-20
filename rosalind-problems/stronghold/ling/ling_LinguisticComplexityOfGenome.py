import math
from BioInfoToolkit.Sequences.SuffixTree import SuffixTreeUkkonen
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/ling/

Given a length n string s formed over an alphabet A of size a, let the "substring count" sub(s) denote the total number of distinct substrings of s. Furthermore, let the "maximum substring count" m(a,n) denote the maximum number of distinct substrings that could appear in a string of length n formed over A.

The linguistic complexity of s (written lc(s)) is equal to sub(s)m(a,n); in other words, lc(s) represents the percentage of observed substrings of s to the total number that are theoretically possible. Note that 0<lc(s)<1, with smaller values of lc(s) indicating that s is more repetitive.

As an example, consider the DNA string (a=4) s=ATTTGGATT. In the following table, we demonstrate that lc(s)=3540=0.875 by considering the number of observed and possible length k substrings of s, which are denoted by subk(s) and m(a,k,n), respectively. (Observe that m(a,n)=∑nk=1m(a,k,n)=40 and sub(s)=∑nk=1subk(s)=35.)

    Given: A DNA string s of length at most 100 kbp.

    Return: The linguistic complexity lc(s).
"""


def lc(a: int, seq: str):
    def m_a_k(a: int, n: int, k: int):
        if k == 1:
            return a
        if k*math.log(a) < math.log(n-k+1):
            return a**k
        return n-k+1

    def sub(seq: str) -> int:
        sTree = SuffixTreeUkkonen(seq+'$')
        uSubCount = sTree.countUniqueSubstrings()
        return uSubCount

    n = len(seq)
    m_a_n = sum(m_a_k(a, n, k) for k in range(1, n+1))
    uSub = sub(seq)

    return uSub/m_a_n



def verify(result: float, solution: float) -> bool:
    correct = abs(result-solution) <= 0.001
    return correct


def solve(seq: str) -> float:
    a = 4
    lc_res = lc(a, seq)
    return lc_res


def load_results(path: str) -> float:
    lines = readTextFile(path)
    lc_res = float(lines[0])
    return lc_res


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq = lines[0]

    result = solve(seq)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ling_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    lc_res = solve(seq)

    print(lc_res)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(lc_res), 'w')

    correct = solve_and_check(path)
    print(correct)
