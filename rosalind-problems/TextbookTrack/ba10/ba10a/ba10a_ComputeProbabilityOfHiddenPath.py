import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.IO.IO_utils import matrix_text_to_dict, trans_matrix_from_dict

"""
https://rosalind.info/problems/ba10a/

Probability of a Hidden Path Problem

    Given: A hidden path π followed by the states States and transition matrix Transition of an HMM (Σ, States, Transition, Emission).

    Return: The probability of this path, Pr(π). You may assume that initial probabilities are equal.
"""

OutputT = float


def compute_path_prob(path: str, states: list[str], trans_dict: dict[tuple[str, str], float]):
    p0 = 1/len(states)
    prob = p0
    for a, b in zip(path, path[1:]):
        pt = trans_dict[(a, b)]
        prob *= pt
    return prob


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = abs(result-solution) / solution <= 0.001
    return correct


def solve(seq: str, state_space: list[str], trans_dict: dict[tuple[str, str], float]) -> OutputT:
    prob_path = compute_path_prob(seq, state_space, trans_dict)
    return prob_path


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    prob = float(lines[0])
    return prob


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq = lines[0]
    state_space = lines[2].split()
    n = len(state_space)
    trans_matrix_header = lines[4].split()
    trans_matrix_lines = [line.split() for line in lines[5:5+n]]

    trans_dict = matrix_text_to_dict(trans_matrix_header, trans_matrix_lines)
    # trans_matrix = trans_matrix_from_dict(trans_dict, state_space)

    result = solve(seq, state_space, trans_dict)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba10a_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    state_space = lines[2].split()
    n = len(state_space)
    trans_matrix_header = lines[4].split()
    trans_matrix_lines = [line.split() for line in lines[5:5+n]]

    trans_dict = matrix_text_to_dict(trans_matrix_header, trans_matrix_lines)
    trans_matrix = trans_matrix_from_dict(trans_dict, state_space)

    prob_path = solve(seq, state_space, trans_dict)

    out = str(prob_path)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

