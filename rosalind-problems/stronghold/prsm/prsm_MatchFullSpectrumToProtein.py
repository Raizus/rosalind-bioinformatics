from collections import Counter
import BioInfoToolkit.Spectrometry as Spectrometry
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMass


"""
https://rosalind.info/problems/prsm/

The complete spectrum of a weighted string s is the multiset S[s] containing the weights of every prefix and suffix of s.

    Given: A positive integer n followed by a collection of n protein strings s1, s2, ..., sn and a multiset R of positive numbers (corresponding to the complete spectrum of some unknown protein string).

    Return: The maximum multiplicity of R⊖S[sk] taken over all strings sk, followed by the string sk for which this maximum multiplicity occurs (you may output any such value if multiple solutions exist).
"""


def verify(result: tuple[int, str], solution: tuple[int, str]) -> bool:
    mult_res = result[0]
    mult_sol = solution[0]

    correct = mult_res == mult_sol

    return correct


def solve(proteins: list[str], spectrum: Counter[float]) -> tuple[int, list[str]]:
    multiplicityDict: dict[str, int] = dict()
    for _, prot in enumerate(proteins):
        Ssk = Counter(Spectrometry.getCompleteMassSpectrum(
            prot, AminoacidMonoisotopicMass))
        resultMs = Spectrometry.minkowskiDiff(spectrum, Ssk)
        maxMultiplicity = Spectrometry.multisetMaxMultiplicity(resultMs)
        multiplicityDict[prot] = maxMultiplicity

    max_mult = max(multiplicityDict.values())
    strings = [key for key, val in multiplicityDict.items()
                   if val == max_mult]
    
    return max_mult, strings


def load_results(path: str) -> tuple[int, str]:
    lines = readTextFile(path)
    max_mult = int(lines[0])
    string = lines[1]

    return max_mult, string


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    n = int(lines[0])
    proteins = lines[1:n+1]
    spectrum = Counter([float(val) for val in lines[n+1:]])

    (mult, strings) = solve(proteins, spectrum)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify((mult, strings[0]), solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_prsm_1.txt'

    lines = readTextFile(path)
    n = int(lines[0])
    protStrs = lines[1:n+1]
    spectrum = Counter([float(val) for val in lines[n+1:]])

    multiplicityDict: dict[str, int] = dict()
    for i, prot in enumerate(protStrs):
        Ssk = Counter(Spectrometry.getCompleteMassSpectrum(
            prot, AminoacidMonoisotopicMass))
        resultMs = Spectrometry.minkowskiDiff(spectrum, Ssk)
        maxMultiplicity = Spectrometry.multisetMaxMultiplicity(resultMs)
        multiplicityDict[prot] = maxMultiplicity

    maxMult = max(multiplicityDict.values())
    maxMultStrs = [key for key, val in multiplicityDict.items()
                   if val == maxMult]
    
    print(maxMult)
    print(maxMultStrs[0])

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(maxMult), 'w')
    writeTextFile(result_path, maxMultStrs[0], 'a')

    # correct = solve_and_check(path)
    # print(correct)
