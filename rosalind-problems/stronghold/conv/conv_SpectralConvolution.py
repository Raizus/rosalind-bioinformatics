from collections import Counter
import BioInfoToolkit.Spectrometry as Spectrometry
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/conv/

A multiset is a generalization of the notion of set to include a collection of objects in which each object may occur more than once (the order in which objects are given is still unimportant). For a multiset S, the multiplicity of an element x is the number of times that x occurs in the set; this multiplicity is denoted S(x). Note that every set is included in the definition of multiset.

The Minkowski sum of multisets S1 and S2 containing real numbers is the new multiset S1⊕S2 formed by taking all possible sums s1+s2 of an element s1 from S1 and an element s2 from S2. The Minkowski sum could be defined more concisely as S1⊕S2=s1+s2:s1∈S1,s2∈S2, The Minkowski difference S1⊖S2 is defined analogously by taking all possible differences s1-s2.

If S1 and S2 represent simplified spectra taken from two peptides, then S1⊖S2 is called the spectral convolution of S1 and S2. In this notation, the shared peaks count is represented by (S2⊖S1)(0), and the value of x for which (S2⊖S1)(x) has the maximal value is the shift value maximizing the number of shared masses of S1 and S2.

    Given: Two multisets of positive real numbers S1 and S2. The size of each multiset is at most 200.

    Return: The largest multiplicity of S1⊖S2, as well as the absolute value of the number x maximizing (S1⊖S2)(x) (you may return any such value if multiple solutions exist).
"""


def verify(result: tuple[int, float], solution: tuple[int, float], 
           ms1: Counter[float], ms2: Counter[float]) -> bool:
    mult_res, x_res = result
    mult_sol, x_sol = solution

    newSet = Spectrometry.minkowskiDiff(ms1, ms2, 8)
    maxCount = max(newSet.values())

    correct = mult_res == mult_sol == maxCount and newSet[x_res] == newSet[x_sol]

    return correct


def solve(ms1: Counter[float], ms2: Counter[float]) -> tuple[int, float]:
    newSet = Spectrometry.minkowskiDiff(ms1, ms2, 8)
    maxCount = max(newSet.values())
    maxKeys = [key for key, val in newSet.items() if val == maxCount]

    return maxCount, maxKeys[0]


def load_results(path: str) -> tuple[int, float]:
    lines = readTextFile(path)
    mult = int(lines[0])
    x = float(lines[1])

    return mult, x


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    ms1 = Counter([float(val) for val in lines[0].split(" ")])
    ms2 = Counter([float(val) for val in lines[1].split(" ")])
    result = solve(ms1, ms2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, ms1, ms2)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_conv_1.txt'

    # If S1 and S2 represent simplified spectra taken from two peptides, then S1 ⊖ S2 is called the spectral convolution of S1 and S2. In this notation, the shared peaks count is represented by (S2⊖S1)(0), and the value of x for which (S2⊖S1)(x) has the maximal value is the shift value maximizing the number of shared masses of S1 and S2.

    lines = readTextFile(path)
    ms1 = Counter([float(val) for val in lines[0].split(" ")])
    ms2 = Counter([float(val) for val in lines[1].split(" ")])
    mult, x = solve(ms1, ms2)

    print(mult)
    print(x)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(mult), 'w')
    writeTextFile(result_path, str(x), 'a')

    correct = solve_and_check(path)
    print(correct)
