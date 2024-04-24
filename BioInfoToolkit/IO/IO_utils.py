import numpy as np

from BioInfoToolkit.IO.IO import writeTextFile

def matrix_text_to_dict(matrix_header: list[str], matrix_lines: list[list[str]]) -> dict[tuple[str, str], float]:
    matrix_dict: dict[tuple[str, str], float] = dict()
    for line in matrix_lines:
        state2 = line[0]
        for i, state in enumerate(matrix_header):
            prob = float(line[i+1])
            matrix_dict[(state2, state)] = prob

    return matrix_dict

def trans_matrix_from_dict(trans_dict: dict[tuple[str, str], float], state_space: list[str]):
    K = len(state_space)
    trans_matrix = np.zeros((K,K), dtype=np.float64)
    for i, state1 in enumerate(state_space):
        for j, state2 in enumerate(state_space):
            trans_matrix[i,j] = trans_dict.get((state1, state2), 0)
    return trans_matrix

def emission_matrix_from_dict(emission_dict: dict[tuple[str, str], float], state_space: list[str], obs_space: list[str]):
    K = len(state_space)
    T = len(obs_space)
    emission_matrix = np.zeros((K,T), dtype=np.float64)
    for i, state1 in enumerate(state_space):
        for j, obs in enumerate(obs_space):
            emission_matrix[i,j] = emission_dict.get((state1, obs), 0)
    return emission_matrix


def print_and_save_matrix(row_header: list[str], col_header: list[str], matrix, path_res: str, precision: int = 3, spacing: str = '\t'):
    header = spacing+spacing.join(row_header)
    writeTextFile(path_res, header, 'a')
    print(header)
    for i, line in enumerate(matrix):
        text = f"{col_header[i]}{spacing}" + \
            spacing.join(str(round(val, precision)) for val in line)
        print(text)
        writeTextFile(path_res, text, 'a')
