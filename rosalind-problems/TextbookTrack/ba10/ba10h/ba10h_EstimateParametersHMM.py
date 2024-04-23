from BioInfoToolkit.HMM import estimateParametersHMM
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.IO.IO_utils import print_and_save_matrix
from testing import float_matrices_match

"""
https://rosalind.info/problems/ba10h/

HMM Parameter Estimation Problem

    Given: A sequence of emitted symbols x = x1 . . . xn in an alphabet ∑ and a path π = π1 . . . πn generated by a k-state HMM with unknown transition and emission probabilities.

    Return: A matrix of transition probabilities Transition and a matrix of emission probabilities Emission that maximize Pr(x,π) over all possible matrices of transition and emission probabilities.
"""

OutputT = tuple[list[list[float]], list[list[float]]]


def verify(result: OutputT, solution: OutputT) -> bool:
    trans_matrix_res, emission_matrix_res = result
    trans_matrix_sol, emission_matrix_sol = solution
    tm = float_matrices_match(trans_matrix_res, trans_matrix_sol)
    em = float_matrices_match(emission_matrix_res, emission_matrix_sol)
    correct = tm and em
    return correct


def solve(observations: str, obs_space: list[str], hidden_path: str, state_space: list[str]) -> OutputT:
    state_space_map: dict[str, int] = {state: i for i, state in enumerate(state_space)}
    obs_space_map: dict[str, int] = {obs: i for i, obs in enumerate(obs_space)}

    hidden_path_list = [c for c in hidden_path]
    trans_matrix, emission_matrix = estimateParametersHMM(
        observations, hidden_path_list, state_space_map, obs_space_map)

    trans_matrix_: list[list[float]] = trans_matrix.tolist()
    emission_matrix_: list[list[float]] = emission_matrix.tolist()

    return trans_matrix_, emission_matrix_


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

    return trans_matrix, emission_matrix


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    observations = lines[0]
    obs_space = lines[2].split()
    hidden_path = lines[4]
    state_space = lines[6].split()

    result = solve(
        observations, obs_space, hidden_path, state_space)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba10h_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    observations = lines[0]
    obs_space = lines[2].split()
    hidden_path = lines[4]
    state_space = lines[6].split()

    trans_matrix, emission_matrix = solve(
        observations, obs_space, hidden_path, state_space)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    spacing = '\t'

    print_and_save_matrix(state_space, state_space,
                          trans_matrix, result_path, precision=4)

    print('-'*8)
    writeTextFile(result_path, '-'*8, 'a')

    print_and_save_matrix(obs_space, state_space,
                          emission_matrix, result_path, precision=4)

    correct = solve_and_check(path)
    print(correct)