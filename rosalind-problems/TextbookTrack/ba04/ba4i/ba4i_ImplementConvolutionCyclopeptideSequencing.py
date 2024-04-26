from BioInfoToolkit.Spectrometry.Spectrometry import convolutionLeaderboardCyclopeptideSequencing, peptideSpectrumScoreTuple
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4i/

Implement ConvolutionCyclopeptideSequencing
    Given: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.

    Return: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
"""


OutputT = str


def verify(result: OutputT, solution: OutputT, spectrum: list[int]) -> bool:
    peptide_res = tuple(int(v) for v in result.split('-'))
    peptide_sol = tuple(int(v) for v in solution.split('-'))

    score_res = peptideSpectrumScoreTuple(
        peptide_res, spectrum, True)
    score_sol = peptideSpectrumScoreTuple(
        peptide_sol, spectrum, True)

    correct = score_res <= score_sol
    return correct


def solve(M: int, N: int, spectrum: list[int]) -> OutputT:
    # mass_dict = AminoacidMonoisotopicMassInteger
    leader_peptide, _ = convolutionLeaderboardCyclopeptideSequencing(
        spectrum, M, N)

    peptide = '-'.join(str(aa) for aa in leader_peptide)

    return peptide


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptide = lines[0]
    return peptide


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    M = int(lines[0])
    N = int(lines[1])
    spectrum = [int(m) for m in lines[2].split()]
    spectrum = sorted(spectrum)

    result = solve(M, N, spectrum)

    solution = load_results(solution_path)

    correct = verify(result, solution, spectrum)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4i_3.txt'

    lines = readTextFile(path)
    M = int(lines[0])
    N = int(lines[1])
    spectrum = [int(m) for m in lines[2].split()]
    spectrum = sorted(spectrum)

    # Note: Highly probable that there's a mistake in rosalind's verifier

    peptide = solve(M, N, spectrum)

    out = str(peptide)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
