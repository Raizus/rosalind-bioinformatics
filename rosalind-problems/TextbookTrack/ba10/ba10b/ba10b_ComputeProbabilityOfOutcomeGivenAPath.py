import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.IO.IO_utils import matrix_text_to_dict, trans_matrix_from_dict

"""
https://rosalind.info/problems/ba10b/

Probability of an Outcome Given a Hidden Path Problem

    Given: A string x, followed by the alphabet Σ from which x was constructed, followed by a hidden path π, followed by the states States and emission matrix Emission of an HMM (Σ, States, Transition, Emission).

    Return: The conditional probability Pr(x|π) that string x will be emitted by the HMM given the hidden path π.
"""

OutputT = float


def compute_conditional_string_prob(x: str, path_pi: str, trans_dict: dict[tuple[str, str], float]):
    prob = 1
    for a, b in zip(path_pi, x):
        pt = trans_dict[(a, b)]
        prob *= pt
    return prob


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = abs(result-solution) / solution <= 0.001
    return correct


def solve(x: str, path_pi: str, trans_dict: dict[tuple[str, str], float]) -> OutputT:
    prob_path = compute_conditional_string_prob(seq, path_pi, trans_dict)
    return prob_path


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    prob = float(lines[0])
    return prob


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seq = lines[0]
    alphabet = lines[2].split()
    n = len(alphabet)
    path_pi = lines[4]
    state_space = lines[6].split()
    trans_matrix_header = lines[8].split()
    trans_matrix_lines = [line.split() for line in lines[9:9+n]]

    trans_dict = matrix_text_to_dict(trans_matrix_header, trans_matrix_lines)
    # trans_matrix = trans_matrix_from_dict(trans_dict, state_space)

    result = solve(seq, path_pi, trans_dict)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba10b_1.txt'

    lines = readTextFile(path)
    seq = lines[0]
    alphabet = lines[2].split()
    n = len(alphabet)
    path_pi = lines[4]
    state_space = lines[6].split()
    trans_matrix_header = lines[8].split()
    trans_matrix_lines = [line.split() for line in lines[9:9+n]]

    trans_dict = matrix_text_to_dict(trans_matrix_header, trans_matrix_lines)
    trans_matrix = trans_matrix_from_dict(trans_dict, state_space)

    prob_path = solve(seq, path_pi, trans_dict)

    out = str(prob_path)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
