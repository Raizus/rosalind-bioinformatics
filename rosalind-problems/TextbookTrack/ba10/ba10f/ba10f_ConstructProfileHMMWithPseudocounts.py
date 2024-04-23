from BioInfoToolkit.HMM import profileHMM
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.IO.IO_utils import print_and_save_matrix
from testing import float_matrices_match

"""
https://rosalind.info/problems/ba10f/

Profile HMM with Pseudocounts Problem

    Given: A threshold θ and a pseudocount σ, followed by an alphabet Σ, followed by a multiple alignment Alignment whose strings are formed from Σ.

    Return: The transition and emission probabilities of the profile HMM HMM(Alignment, θ, σ).
"""

OutputT = tuple[list[str], list[list[float]], list[list[float]]]


def verify(result: OutputT, solution: OutputT) -> bool:
    state_space_res, trans_matrix_res, emission_matrix_res = result
    state_space_sol, trans_matrix_sol, emission_matrix_sol = solution
    tm = float_matrices_match(trans_matrix_res, trans_matrix_sol)
    em = float_matrices_match(emission_matrix_res, emission_matrix_sol)
    correct = tm and em and state_space_res == state_space_sol
    return correct


def solve(alignment: list[str], obs_space: list[str], threshold: float, pseudocounts: float) -> OutputT:
    trans_matrix, emission_matrix, state_space = profileHMM(
        alignment, obs_space, threshold, pseudocounts)

    trans_matrix_: list[list[float]] = trans_matrix.tolist()
    emission_matrix_: list[list[float]] = emission_matrix.tolist()

    return state_space, trans_matrix_, emission_matrix_


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    state_space = lines[0].split()
    n = len(state_space)

    trans_matrix: list[list[float]] = []
    for line in lines[1:1+n]:
        vals = line.split()
        vals = [float(v) for v in vals[1:]]
        trans_matrix.append(vals)

    emission_matrix: list[list[float]] = []
    for line in lines[1+n+2:1+n+2+n]:
        vals = line.split()
        vals = [float(v) for v in vals[1:]]
        emission_matrix.append(vals)

    return state_space, trans_matrix, emission_matrix


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    threshold, pseudocounts = [float(a) for a in lines[0].split()]
    alphabet = lines[2].split()
    alignment = lines[4:]

    result = solve(alignment, alphabet, threshold, pseudocounts)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba10f_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    threshold, pseudocounts = [float(a) for a in lines[0].split()]
    alphabet = lines[2].split()
    alignment = lines[4:]

    state_space, trans_matrix, emission_matrix = solve(
        alignment, alphabet, threshold, pseudocounts)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    spacing = '\t'

    print_and_save_matrix(state_space, state_space,
                          trans_matrix, result_path, precision=4)

    print('-'*8)
    writeTextFile(result_path, '-'*8, 'a')

    print_and_save_matrix(alphabet, state_space,
                          emission_matrix, result_path, precision=4)

    correct = solve_and_check(path)
    print(correct)
