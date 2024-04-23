import numpy as np
import numpy.typing as npt
from BioInfoToolkit.HMM import sequenceAlignmentWithProfileHMM
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


# https://www.youtube.com/watch?v=HbA0odlLuZs
"""
https://rosalind.info/problems/ba10g/

Sequence Alignment with Profile HMM Problem

    Given: A string Text, a multiple alignment Alignment, a threshold θ, and a pseudocount σ.

    Return: An optimal hidden path emitting Text in HMM(Alignment,θ,σ).
"""


def alignmentWithProfileHMMEvaluate(observations: str, hidden_states: list[str],
                                    state_space: list[str],
                                    trans_matrix: npt.NDArray[np.float64],
                                    emission_matrix: npt.NDArray[np.float64]):
    observation_map = {o: i for i, o in enumerate(alphabet)}
    state_space_map = {s: i for i, s in enumerate(state_space)}

    log_p = 0
    for i, state in enumerate(hidden_states[1:-1], 1):
        si = state_space_map[state]
        psi = state_space_map[hidden_states[i-1]]
        tp = trans_matrix[psi, si]

        obs = observations[i-1]
        oi = observation_map.get(obs)
        st = state[0]
        ep = 1 if obs == '-' and st == 'D' else emission_matrix[si,
                                                                oi] if oi is not None else 0

        log_p += np.log(tp) + np.log(ep)

    si = state_space_map['E']
    psi = state_space_map[hidden_states[-2]]
    log_p += np.log(trans_matrix[psi, si])
    return log_p


OutputT = list[str]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(alignment: list[str], obs_space: list[str], observations: str, threshold: float, pseudocounts: float) -> OutputT:
    hidden_states, _, _ = sequenceAlignmentWithProfileHMM(
        alignment, obs_space, observations, threshold, pseudocounts)

    return hidden_states[1:-1]


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    hidden_states = lines[0].split()

    return hidden_states


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    observations = lines[0]
    threshold, pseudocounts = [float(a) for a in lines[2].split()]
    alphabet = lines[4].split()
    alignment = lines[6:]

    result = solve(alignment, alphabet, observations, threshold, pseudocounts)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba10g_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]

    observations = lines[0]
    threshold, pseudocounts = [float(a) for a in lines[2].split()]
    alphabet = lines[4].split()
    alignment = lines[6:]

    hidden_states = solve(
        alignment, alphabet, observations, threshold, pseudocounts)

    out = ' '.join(hidden_states)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
