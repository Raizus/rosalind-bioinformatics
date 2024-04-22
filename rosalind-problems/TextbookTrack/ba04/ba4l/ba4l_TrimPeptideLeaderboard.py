from collections import Counter
from BioInfoToolkit.IO.IO import readTextFile
from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, peptideSpectrumScore
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4l/

The Trim algorithm, shown below, sorts all peptides in Leaderboard according to their scores, resulting in a sorted Leaderboard. Trim> then retains the top N scoring peptides including ties, and removes all other peptides from Leaderboard.

    Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
        for j ← 1 to |Leaderboard|
            Peptide ← j-th peptide in Leaderboard
            LinearScores(j) ← LinearScore(Peptide, Spectrum)
        sort Leaderboard according to the decreasing order of scores in LinearScores
        sort LinearScores in decreasing order
        for j ← N + 1 to |Leaderboard|
            if LinearScores(j) < LinearScores(N)
                remove all peptides starting from the j-th peptide from Leaderboard
            return Leaderboard
        return Leaderboard

Trim Problem
    Trim a leaderboard of peptides.

        Given: A leaderboard of linear peptides Leaderboard, a linear spectrum Spectrum, and an integer N.

        Return: The top N peptides from Leaderboard scored against Spectrum. Remember to use LinearScore.
"""


def trim_leaderboard(leaderboard: list[str], spectrum: list[int], N: int, massDict: dict[str, int]):
    if len(leaderboard) <= N:
        return leaderboard

    scores: list[tuple[str, int]] = []
    for peptide in leaderboard:
        score = peptideSpectrumScore(peptide, spectrum, massDict, False)
        scores.append((peptide, score))
    scores.sort(key=lambda x: x[1], reverse=True)

    scoreN = scores[N-1][1]
    leaderboard = [pep for pep, score in scores if score >= scoreN]
    return leaderboard

OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = Counter(result) == Counter(solution)
    return correct


def solve(leaderboard: list[str], spectrum: list[int], N: int) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger
    leaderboard = trim_leaderboard(leaderboard, spectrum, N, massDict)

    return leaderboard


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    leaderboard = lines[0].split()
    return leaderboard


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    leaderboard = lines[0].split()
    spectrum = [int(v) for v in lines[1].split(' ')]
    N = int(lines[2])

    result = solve(leaderboard, spectrum, N)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4l_1.txt'

    lines = readTextFile(path)
    leaderboard = lines[0].split()
    spectrum = [int(v) for v in lines[1].split(' ')]
    N = int(lines[2])

    leaderboard2 = solve(leaderboard, spectrum, N)

    out = ' '.join(leaderboard2)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
