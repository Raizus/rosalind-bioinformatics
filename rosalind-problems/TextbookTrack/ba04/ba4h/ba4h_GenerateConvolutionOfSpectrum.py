from collections import Counter
from BioInfoToolkit.Spectrometry.Spectrometry import spectrumConvolution
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import math

"""
https://rosalind.info/problems/ba4h/

We define the convolution of a cyclic spectrum by taking all positive differences of masses in the spectrum. Figure 1 shows the convolution of the theoretical spectrum of "NQEL".

As predicted, some of the values in Figure 2.12 appear more frequently than others. For example, 113 (the mass of "L") appears eight times in the convolution of the theoretical spectrum of "NQEL"; we say that 113 has multiplicity 8. Six of the eight occurrences of 113 correspond to subpeptide pairs differing in an "L": "L" and ""; "LN" and "N"; "EL" and "E"; "LNQ" and "NQ"; "QEL" and "QE"; "NQEL" and "NQE".

Spectral Convolution Problem
    Compute the convolution of a spectrum.

        Given: A collection of integers Spectrum.

        Return: The list of elements in the convolution of Spectrum in decreasing order of their multiplicities. If an element has multiplicity k, it should appear exactly k times.
"""


def sortByMultiplicity(spectrum_conv: Counter[int]):
    sorted_counter = sorted(
        (item for item in spectrum_conv.items()), key=lambda x: x[1], reverse=True)
    sorted_spectrum: list[int] = []
    for mass, count in sorted_counter:
        sorted_spectrum.extend([mass]*count)
    return sorted_spectrum


OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    ms_res = Counter(result)
    if ms_res != Counter(solution):
        return False
    
    # check if result is correctly sorted
    last = math.inf
    for m in result:
        mult = ms_res[m]
        if mult > last:
            return False

    return True


def solve(spectrum: list[int]) -> OutputT:
    conv = spectrumConvolution(spectrum)
    conv = sortByMultiplicity(conv)

    return conv


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    conv = [int(v) for v in lines[0].split()]
    return conv


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    spectrum = [int(m) for m in lines[0].split()]

    result = solve(spectrum)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4h_1.txt'

    lines = readTextFile(path)
    spectrum = [int(m) for m in lines[0].split()]

    conv = solve(spectrum)

    out = ' '.join(str(v) for v in conv)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
