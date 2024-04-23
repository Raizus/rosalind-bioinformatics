from BioInfoToolkit.HMM import forward_backward, viterbi_learning
from BioInfoToolkit.IO.IO import readTextFile, writeTextFile
import numpy as np
from IO_utils import matrix_text_to_dict, trans_matrix_from_dict, emission_matrix_from_dict, print_and_save_trans_matrix, print_and_save_emission_matrix


if __name__ == "__main__":
    path = 'rosalind-problems/TextbookTrack/test_data/test_ba10j_3.txt'
    # path = 'rosalind-problems/TextbookTrack/test_data/rosalind_ba10j.txt'

    lines = readTextFile(path)

    observations = lines[0]
    obs_space = lines[2].split()
    state_space = lines[4].split()

    K = len(state_space)
    T = len(obs_space)

    trans_matrix_header = lines[6].split()
    trans_matrix_lines = [line.split() for line in lines[7:7+K]]

    idx = 7+K+1
    emission_matrix_header = lines[idx].split()
    emission_matrix_lines = [line.split() for line in lines[idx+1:]]

    trans_dict = matrix_text_to_dict(trans_matrix_header, trans_matrix_lines)
    emission_dict = matrix_text_to_dict(
        emission_matrix_header, emission_matrix_lines)
    trans_matrix = trans_matrix_from_dict(trans_dict, state_space)
    emission_matrix = emission_matrix_from_dict(
        emission_dict, state_space, obs_space)

    initial_probs = np.array([1/len(state_space) for _ in state_space])
    
    forward_backward(state_space, obs_space, observations, initial_probs, trans_matrix, emission_matrix)
