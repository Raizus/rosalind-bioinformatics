from BioInfoToolkit.HMM import HMM_forward_algorithm
import numpy as np
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.IO.IO_utils import emission_matrix_from_dict, matrix_text_to_dict, trans_matrix_from_dict
import numpy.typing as npt

"""
https://rosalind.info/problems/ba10d/

Outcome Likelihood Problem

    Given: A string x, followed by the alphabet Σ from which x was constructed, followed by the states States, transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).

    Return: The probability Pr(x) that the HMM emits x. 
"""

OutputT = float


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = abs(result-solution) / solution <= 0.001
    return correct


def solve(obs_space: list[str], state_space: list[str],
          observations: str, initial_probs: npt.NDArray,
          trans_matrix: npt.NDArray, emission_matrix: npt.NDArray) -> OutputT:
    prob = HMM_forward_algorithm(
        obs_space, state_space, observations, initial_probs, trans_matrix, emission_matrix)

    return prob


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    prob = float(lines[0])
    return prob


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)

    observations = lines[0]

    obs_space = lines[2].split()
    N = len(obs_space)
    state_space = lines[4].split()
    K = len(state_space)

    trans_matrix_header = lines[6].split()
    trans_matrix_lines = [line.split() for line in lines[7:7+K]]

    idx = 7+K+1
    emission_matrix_header = lines[idx].split()
    emission_matrix_lines = [line.split() for line in lines[idx+1:idx+1+N]]

    initial_probs = np.array([1/len(state_space) for _ in state_space])

    trans_dict = matrix_text_to_dict(trans_matrix_header, trans_matrix_lines)
    emission_dict = matrix_text_to_dict(
        emission_matrix_header, emission_matrix_lines)
    trans_matrix = trans_matrix_from_dict(trans_dict, state_space)
    emission_matrix = emission_matrix_from_dict(
        emission_dict, state_space, obs_space)

    result = solve(obs_space, state_space, observations,
                   initial_probs, trans_matrix, emission_matrix)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba10d_1.txt'

    lines = readTextFile(path)

    observations = lines[0]

    obs_space = lines[2].split()
    N = len(obs_space)
    state_space = lines[4].split()
    K = len(state_space)

    trans_matrix_header = lines[6].split()
    trans_matrix_lines = [line.split() for line in lines[7:7+K]]

    idx = 7+K+1
    emission_matrix_header = lines[idx].split()
    emission_matrix_lines = [line.split() for line in lines[idx+1:idx+1+N]]

    initial_probs = np.array([1/len(state_space) for _ in state_space])

    trans_dict = matrix_text_to_dict(trans_matrix_header, trans_matrix_lines)
    emission_dict = matrix_text_to_dict(
        emission_matrix_header, emission_matrix_lines)
    trans_matrix = trans_matrix_from_dict(trans_dict, state_space)
    emission_matrix = emission_matrix_from_dict(
        emission_dict, state_space, obs_space)

    prob = solve(obs_space, state_space, observations,
                      initial_probs, trans_matrix, emission_matrix)

    out = str(prob)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
