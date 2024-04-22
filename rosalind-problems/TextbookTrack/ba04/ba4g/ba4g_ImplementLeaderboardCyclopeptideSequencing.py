from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, leaderboardCyclopeptideSequencing, peptideSpectrumScoreTuple
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba4g/

We have thus far worked with theoretical spectra of cyclic peptides, in which the mass of every subpeptide is given. This inflexibility presents a practical barrier, since mass spectrometers generate spectra that are far from ideal — they are characterized by having both false masses and missing masses. A false mass is present in the experimental spectrum but absent from the theoretical spectrum; a missing mass is present in the theoretical spectrum but absent from the experimental spectrum (see Figure 1).

To generalize “Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum” to handle “noisy” spectra having false and missing masses, we need to relax the requirement that a candidate peptide's theoretical spectrum must match the experimental spectrum exactly, and instead incorporate a scoring function that will select the peptide whose theoretical spectrum matches the given experimental spectrum the most closely. Given a cyclic peptide Peptide and a spectrum Spectrum, we define Score(Peptide, Spectrum) as the number of masses shared between Cyclospectrum(Peptide) and Spectrum. Recalling Figure 1, if

    Spectrum = {0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484},

then Score("NQEL", Spectrum) = 11.

To limit the number of candidate peptides under consideration, we will use a Leaderboard, which holds the N highest scoring candidate peptides for further extension. At each step, we will expand all candidate peptides found in Leaderboard by adding every possible amino acid to the end. Then, we will eliminate those peptides whose newly calculated scores are not high enough to keep them on the Leaderboard. This idea is similar to the notion of a “cut” in a golf tournament; after the cut, only the top N golfers are allowed to play in the next round, since they are the only players who have a reasonable chance of winning.

To be fair, a cut should include anyone who is tied with the Nth-place competitor. Thus, Leaderboard should be trimmed down to the “N highest-scoring peptides including ties”, which may include more than N peptides. Given a list of peptides Leaderboard, a spectrum Spectrum, and an integer N, Cut(Leaderboard, Spectrum, N) returns the top N highest-scoring peptides in <Leaderboard (including ties) with respect to Spectrum. We now introduce LEADERBOARDCYCLOPEPTIDESEQUENCING. In what follows, the 0-peptide is the peptide "" containing no amino acids.

    LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)
        Leaderboard ← {0-peptide}
        LeaderPeptide ← 0-peptide
        while Leaderboard is non-empty
            Leaderboard ← Expand(Leaderboard)
            for each Peptide in Leaderboard
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                        LeaderPeptide ← Peptide
                else if Mass(Peptide) > ParentMass(Spectrum)
                    remove Peptide from Leaderboard
            Leaderboard ← Cut(Leaderboard, Spectrum, N)
        output LeaderPeptide

Implement LeaderboardCyclopeptideSequencing

    Given: An integer N and a collection of integers Spectrum.

    Return: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).
"""


OutputT = str

def verify(result: OutputT, solution: OutputT, spectrum: list[int]) -> bool:
    peptide_res = tuple(int(v) for v in result.split('-'))
    peptide_sol = tuple(int(v) for v in solution.split('-'))

    score_res = peptideSpectrumScoreTuple(
        peptide_res, spectrum, False)
    score_sol = peptideSpectrumScoreTuple(
        peptide_sol, spectrum, False)

    correct = score_res == score_sol
    return correct


def solve(N: int, spectrum: list[int]) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger
    base_masses = list(set(m for m in massDict.values()))
    leaderPeptide, _ = leaderboardCyclopeptideSequencing(spectrum, N, base_masses)
    # print("Parent mass: ", spectrum[-1])
    # print("\tLeader: ", leaderPeptide)
    # print("\tScore: ", leaderScore)
    # print("\tLeader mass: ", sum(leaderPeptide))
    # print('-'.join(str(aa) for aa in leaderPeptide))

    peptide = '-'.join(str(aa) for aa in leaderPeptide)

    return peptide


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    peptide = lines[0]
    return peptide


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    N = int(lines[0])
    spectrum = [int(m) for m in lines[1].split()]

    result = solve(N, spectrum)

    solution = load_results(solution_path)

    correct = verify(result, solution, spectrum)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba4g_1.txt'

    lines = readTextFile(path)
    N = int(lines[0])
    spectrum = [int(m) for m in lines[1].split()]

    score = solve(N, spectrum)

    out = str(score)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
