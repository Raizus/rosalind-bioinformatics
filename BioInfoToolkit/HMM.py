import numpy as np
import numpy.typing as npt
import re

def HMM_forward_algorithm(obs_space: list[str], state_space: list[str],
                          observations: str, initial_probs: npt.NDArray,
                          trans_matrix: npt.NDArray, emission_matrix: npt.NDArray):
    # prior
    obs_space_dict = {obs: i for i, obs in enumerate(obs_space)}
    obs_array = [obs_space_dict[obs] for obs in observations]
    obs_array = np.array(obs_array)
    alpha_i = initial_probs * emission_matrix[:, obs_array[0]]

    for i, obs in enumerate(obs_array[1:], 1):
        alpha_i = emission_matrix[:, obs] * \
            (np.transpose(trans_matrix) @ alpha_i)

    alpha_T = np.sum(alpha_i)

    return alpha_T


def viterbi_best_path(state_space: list[str], T1: npt.NDArray, T2: npt.NDArray) -> list[str]:
    T = T1.shape[1]

    k = np.argmax(T1[:, T-1])
    best_path: list[str] = list()
    for j in range(T-1, -1, -1):
        best_path.append(state_space[k])
        k = T2[k, j]
    best_path.reverse()

    return best_path


def viterbi_compute_matrices(obs_space: list[str], state_space: list[str],
                             observations: str, initial_probs: npt.NDArray,
                             trans_matrix: npt.NDArray, emission_matrix: npt.NDArray):
    N = len(state_space)
    K = len(observations)

    # Each element of T1[i][j] of T1 stores the probability of the most likely path so far X^ = ( x^_1 , x^_2 , … , x^_j ) with x^_j = s_i that generates Y = (y_1 , y_2 , … , y_j ) observations
    T1 = np.zeros((N, K), dtype=np.float64)
    T2 = np.zeros((N, K), dtype=np.int64)

    j = obs_space.index(observations[0])
    T1[:, 0] = initial_probs * emission_matrix[:, j]

    j = 1
    a = 0

    for j, obs in enumerate(observations[1:], 1):
        for i, state in enumerate(state_space):
            obs_idx = obs_space.index(obs)
            k = np.argmax(T1[:, j-1] * trans_matrix[:, i]
                          * emission_matrix[i, obs_idx])

            T2[i, j] = k
            T1[i, j] = T1[k, j-1] * trans_matrix[k, i] * \
                emission_matrix[i, obs_idx]

    return T1, T2


def viterbi_algorithm(obs_space: list[str], state_space: list[str],
                             observations: str, initial_probs: npt.NDArray,
                             trans_matrix: npt.NDArray, emission_matrix: npt.NDArray):
    T1, T2 = viterbi_compute_matrices(obs_space, state_space, observations, initial_probs, trans_matrix, emission_matrix)
    best_path = viterbi_best_path(state_space, T1, T2)
    return best_path


# https://www.youtube.com/watch?v=vO_6xfLwGao
def profileHMM(alignment: list[str], alphabet: list[str],
               threshold: float, pseudocount: float | None = None):
    """Builds a profile HMM from the multiple alignment of sequences

    Args:
        alignment (list[str]): the multiple alignment
        alphabet (list[str]): The alphabet corresponding to possible observations, in this case corresponding to aminoacids or nucleotides.
        threshold (float): 
        pseudocount (float | None, optional): Add pseudocounts to transition and transmission matrix. Should be a float [0, 1[, or None. Defaults to None.

    Returns:
        trans_matrix, emission_matrix, state_space (tuple[NDArray[float64], NDArray[float64], list[str]]): returns a tuple with the transition matrix, the emission matrix and the state space
    """
    n = len(alignment[0])
    nseqs = len(alignment)

    # determine state transitions from alignment
    match_count = 0
    state_transitions: list[list[tuple[str, str]]] = [
        [('S', '')] for _ in range(nseqs)]
    for i in range(n):
        chars = [seq[i] for seq in alignment]
        insertion_count = chars.count('-')
        insertion_ratio = insertion_count/nseqs
        if insertion_ratio < threshold:
            match_count += 1
            for j, char in enumerate(chars):
                if char != '-':
                    state_transitions[j].append((f"M{match_count}", char))
                else:
                    state_transitions[j].append((f"D{match_count}", char))
        else:
            for j, char in enumerate(chars):
                if char != '-':
                    state_transitions[j].append((f"I{match_count}", char))

    for state_transition in state_transitions:
        state_transition.append(('E', ''))

    mid_states = [
        f"{t}{i+1}" for i in range(match_count) for t in ['M', 'D', 'I']]
    state_space = ['S', 'I0'] + mid_states + ['E']
    state_space_map: dict[str, int] = {
        state: i for i, state in enumerate(state_space)}
    obs_space_map: dict[str, int] = {obs: i for i, obs in enumerate(alphabet)}

    N = len(state_space)
    K = len(alphabet)

    trans_matrix = np.zeros((N, N), dtype=np.float64)
    emission_matrix = np.zeros((N, K), dtype=np.float64)

    for state_sequence in state_transitions:
        for (s1, ch1), (s2, ch2) in zip(state_sequence, state_sequence[1:]):
            transition_idx = (
                state_space_map[s1], state_space_map[s2])
            trans_matrix[transition_idx] += 1
            if ch2 in obs_space_map:
                emission_idx = (state_space_map[s2], obs_space_map[ch2])
                emission_matrix[emission_idx] += 1

    # normalize
    for i, v in enumerate(np.sum(trans_matrix, 1)):
        if v != 0:
            trans_matrix[i, :] /= v

    for i, v in enumerate(np.sum(emission_matrix, 1)):
        if v != 0:
            emission_matrix[i, :] /= v

    if pseudocount is None:
        return trans_matrix, emission_matrix, state_space
    
    # add pseudocounts and renormalize probabilities
    allowed_start_end_trans: list[list[tuple[int, int]]] = [
        [(0, 1), (0, 2), (0, 3)],
        [(1, 1), (1, 2), (1, 3)],
        [(N-4, N-2), (N-3, N-2), (N-2, N-2)],
        [(N-4, N-1), (N-3, N-1), (N-2, N-1)]]

    allowed_mid_trans = ([(3*m1+2+i, (3*m1+4+j)) for j in range(3)]
                         for m1 in range(match_count-1) for i in range(3))

    for line_idxs in allowed_start_end_trans:
        for idx in line_idxs:
            trans_matrix[idx] += pseudocount

    for line_idxs in allowed_mid_trans:
        for idx in line_idxs:
            trans_matrix[idx] += pseudocount

    # renormalise
    for i, v in enumerate(np.sum(trans_matrix, 1)):
        if v != 0:
            trans_matrix[i, :] /= v

    emission_state = [
        1+3*m+i for m in range(match_count) for i in range(2)] + [N-2]

    for state in emission_state:
        total = np.sum(emission_matrix[state, :])
        if total == 0:
            emission_matrix[state, :] = 1/K
        else:
            emission_matrix[state, :] += pseudocount

    # renormalise
    for i, v in enumerate(np.sum(emission_matrix, 1)):
        if v != 0:
            emission_matrix[i, :] /= v

    return trans_matrix, emission_matrix, state_space


def alignedSequenceFromHiddenStates(hidden_states: list[str], observations: str):
    aligned_seq = ''
    i = 0
    for state in hidden_states[1:-1]:
        if state[0] in 'IM':
            aligned_seq += observations[i]
            i += 1
        else:  # state[0] == 'D'
            aligned_seq += '-'
    return aligned_seq


def sequenceAlignmentWithProfileHMM(alignment: list[str], alphabet: list[str], observations: str,
                                    threshold: float, pseudocounts: float | None = None):
    trans_matrix, emission_matrix, state_space = profileHMM(
        alignment, alphabet, threshold, pseudocounts)

    # columns on the viterbi Graph:
    columns = [s for s in state_space if re.match(r'S|M[\d+]', s)]
    M = len(columns)
    N = len(observations) + 1
    observation_map = {o: i for i, o in enumerate(alphabet)}
    state_space_map = {s: i for i, s in enumerate(state_space)}

    SM = np.full((N, M), -np.Inf, dtype=np.float64)
    SD = np.full((N, M), -np.Inf, dtype=np.float64)
    SI = np.full((N, M), -np.Inf, dtype=np.float64)

    SD[0, 0] = 0
    SM[0, 0] = 0
    SI[0, 0] = 0

    def weight(s_i: str, s_f: str, i: int, j: int):
        sii = state_space_map.get(s_i)
        sfi = state_space_map.get(s_f)

        sit = s_i[0]
        sft = s_f[0]

        if sii is None or sfi is None:
            tm = 0
        elif i == 1 and sft in 'MI' and sit != 'D':
            tm = 0
        elif j == 1 and sft in 'MD' and sit != 'I':
            tm = 0
        else:
            tm = trans_matrix[sii, sfi]

        return tm

    for j in range(1, M):
        SD[0, j] = SD[0, j-1] + np.log(trans_matrix[3*(j-1), 3*j])

    obs = observation_map[observations[0]] if 0 < len(observations) else None
    e_prob = emission_matrix[1, obs] if obs is not None else 0
    SI[1, 0] = SI[0, 0] + np.log(trans_matrix[0, 1]) + np.log(e_prob)
    for i in range(2, N):
        obs = observation_map[observations[i-1]] if i - \
            1 < len(observations) else None
        e_prob = emission_matrix[1, obs] if obs is not None else 0
        SI[i, 0] = SI[i-1, 0] + np.log(trans_matrix[1, 1]) + np.log(e_prob)

    obs = observation_map[observations[0]] if 0 < len(observations) else None
    e_prob = emission_matrix[2, obs] if obs is not None else 0
    SM[1, 1] = SM[0, 0] + np.log(trans_matrix[0, 2]) + np.log(e_prob)
    e_prob = emission_matrix[4, obs] if obs is not None else 0
    SI[1, 1] = SD[0, 1] + np.log(trans_matrix[3, 4]) + np.log(e_prob)
    SD[1, 1] = SI[1, 0] + np.log(trans_matrix[1, 3])

    for i in range(1, N):
        for j in range(1, M):
            if i == 1 and j == 1:
                continue

            obs = observation_map[observations[i-1]
                                  ] if i-1 < len(observations) else None

            Ij = state_space_map[f"I{j}"]
            e_prob = emission_matrix[Ij, obs] if obs is not None else 0
            # I(j) -> I(j) I(j) -> O(i-1)
            a = SI[i-1, j] + \
                np.log(weight(f"I{j}", f"I{j}", i, j)) + np.log(e_prob)
            # M(j) -> I(j)   I(j) -> O(i-1)
            b = SM[i-1, j] + \
                np.log(weight(f"M{j}", f"I{j}", i, j)) + np.log(e_prob)
            # D(j) -> I(j) I(j) -> O(i-1)
            c = SD[i-1, j] + \
                np.log(weight(f"D{j}", f"I{j}", i, j)) + np.log(e_prob)
            max_val = max(a, b, c)
            SI[i, j] = max_val

            Mj = state_space_map[f"M{j}"]
            e_prob = emission_matrix[Mj, obs] if obs is not None else 0
            # I(j-1) -> M(j) M(j) -> O(i-1)
            a = SI[i-1, j-1] + \
                np.log(weight(f"I{j-1}", f"M{j}", i, j)) + np.log(e_prob)
            # M(j-1) -> M(j)   M(j) -> O(i-1)
            b = SM[i-1, j-1] + \
                np.log(weight(f"M{j-1}", f"M{j}", i, j)) + np.log(e_prob)
            # D(j-1) -> M(j) M(j) -> O(i-1)
            c = SD[i-1, j-1] + \
                np.log(weight(f"D{j-1}", f"M{j}", i, j)) + np.log(e_prob)
            max_val = max(a, b, c)
            SM[i, j] = max_val

            # I(j-1) -> D(j)
            a = SI[i, j-1] + np.log(weight(f"I{j-1}", f"D{j}", i, j))
            # M(j-1) -> D(j)
            b = SM[i, j-1] + np.log(weight(f"M{j-1}", f"D{j}", i, j))
            # D(j-1) -> D(j)
            c = SD[i, j-1] + np.log(weight(f"D{j-1}", f"D{j}", i, j))

            max_val = max(a, b, c)
            SD[i, j] = max_val

    E = state_space_map['E']
    Mj = state_space_map[f"M{M-1}"]
    Ij = state_space_map[f"I{M-1}"]
    Dj = state_space_map[f"D{M-1}"]

    a = SM[N-1, M-1] + np.log(trans_matrix[Mj, E])
    b = SD[N-1, M-1] + np.log(trans_matrix[Dj, E])
    c = SI[N-1, M-1] + np.log(trans_matrix[Ij, E])

    last_states = [f"M{M-1}", f"D{M-1}", f"I{M-1}"]
    log_p = max(a, b, c)

    idxs = [idx for idx, v in enumerate([a, b, c]) if log_p == v]

    # backtrack
    prev_state = last_states[idxs[0]]
    hidden_states: list[str] = [prev_state, 'E']
    i = N - 1
    j = M - 1

    while i > 0 or j > 0:
        obs = observation_map[observations[i-1]] if 0 <= i-1 < len(observations) else None
        if i == 0 and j == 1:
            prev_state = 'S'
            hidden_states = [prev_state] + hidden_states
            j = j - 1
        elif i == 1 and j == 0:
            prev_state = 'S'
            hidden_states = [prev_state] + hidden_states
            i = i - 1
        elif i == 0:
            prev_state = f"D{j-1}"
            hidden_states = [prev_state] + hidden_states
            j = j - 1
        elif j == 0:
            prev_state = f"I{j}"
            hidden_states = [prev_state] + hidden_states
            i = i - 1
        elif i == 1 and j == 1:
            if prev_state[0] == 'M':
                prev_state = 'S'
                hidden_states = [prev_state] + hidden_states
                i, j = i - 1, j - 1
            elif prev_state[0] == 'D':
                j = j - 1
                prev_state = 'I0'
                hidden_states = [prev_state] + hidden_states
            else:
                i = i - 1
                prev_state = 'D1'
                hidden_states = [prev_state] + hidden_states
        elif prev_state[0] == 'M':
            Mj = state_space_map[f"M{j}"]
            e_prob = emission_matrix[Mj, obs] if obs is not None else 0
            # M(j-1) -> M(j)   M(j) -> O(i-1)
            a = SM[i-1, j-1] + \
                np.log(weight(f"M{j-1}", f"M{j}", i, j)) + np.log(e_prob)
            # D(j-1) -> M(j) M(j) -> O(i-1)
            b = SD[i-1, j-1] + \
                np.log(weight(f"D{j-1}", f"M{j}", i, j)) + np.log(e_prob)
            # I(j-1) -> M(j) M(j) -> O(i-1)
            c = SI[i-1, j-1] + \
                np.log(weight(f"I{j-1}", f"M{j}", i, j)) + np.log(e_prob)

            last_states = [f"M{j-1}", f"D{j-1}", f"I{j-1}"]
            idxs = [idx for idx, v in enumerate([a, b, c]) if SM[i, j] == v]
            prev_state = last_states[idxs[0]]
            hidden_states = [prev_state] + hidden_states
            i, j = i - 1, j - 1
        
        elif prev_state[0] == 'D':
            # M(j-1) -> D(j)
            a = SM[i, j-1] + np.log(weight(f"M{j-1}", f"D{j}", i, j))
            # D(j-1) -> D(j)
            b = SD[i, j-1] + np.log(weight(f"D{j-1}", f"D{j}", i, j))
            # I(j-1) -> D(j)
            c = SI[i, j-1] + np.log(weight(f"I{j-1}", f"D{j}", i, j))

            last_states = [f"M{j-1}", f"D{j-1}", f"I{j-1}"]
            idxs = [idx for idx, v in enumerate([a, b, c]) if SD[i, j] == v]
            prev_state = last_states[idxs[0]]
            hidden_states = [prev_state] + hidden_states
            j = j - 1
        else: # prev_state[0] == 'I'
            Ij = state_space_map[f"I{j}"]
            e_prob = emission_matrix[Ij, obs] if obs is not None else 0
            a = SM[i-1, j] + \
                np.log(weight(f"M{j}", f"I{j}", i, j)) + np.log(e_prob)
            b = SD[i-1, j] + \
                np.log(weight(f"D{j}", f"I{j}", i, j)) + np.log(e_prob)
            c = SI[i-1, j] + \
                np.log(weight(f"I{j}", f"I{j}", i, j)) + np.log(e_prob)

            last_states = [f"M{j}", f"D{j}", f"I{j}"]
            idxs = [idx for idx, v in enumerate([a, b, c]) if SI[i, j] == v]
            prev_state = last_states[idxs[0]]
            hidden_states = [prev_state] + hidden_states
            i = i - 1

    # rebuild aligned sequence from hidden states
    aligned_seq = alignedSequenceFromHiddenStates(hidden_states, observations)

    return hidden_states, aligned_seq, log_p


def estimateParametersHMM(observations: str, hidden_path: list[str], state_space_map: dict[str, int], obs_space_map: dict[str, int]):
    N = len(state_space_map)
    K = len(obs_space_map)
    trans_matrix = np.zeros((N, N), dtype=np.float64)
    emission_matrix = np.zeros((N, K), dtype=np.float64)

    for s1, s2 in zip(hidden_path, hidden_path[1:]):
        trans_idx = (state_space_map[s1], state_space_map[s2])
        trans_matrix[trans_idx] += 1

    for i, v in enumerate(np.sum(trans_matrix, 1)):
        if v != 0:
            trans_matrix[i, :] /= v
        else:
            trans_matrix[i, :] = 1/N

    for state, obs in zip(hidden_path, observations):
        emission_idx = (state_space_map[state], obs_space_map[obs])
        emission_matrix[emission_idx] += 1

    for i, v in enumerate(np.sum(emission_matrix, 1)):
        if v != 0:
            emission_matrix[i, :] /= v
        else:
            emission_matrix[i, :] = 1/K

    return trans_matrix, emission_matrix


def viterbi_learning(state_space: list[str], obs_space: list[str], observations: str, 
                     trans_matrix: npt.NDArray, emission_matrix: npt.NDArray, iters: int):
    state_space_map: dict[str, int] = {state: i for i, state in enumerate(state_space)}
    obs_space_map: dict[str, int] = {obs: i for i, obs in enumerate(obs_space)}

    initial_probs = np.array([1/len(state_space) for _ in state_space])

    for i in range(iters):
        hidden_path = viterbi_algorithm(obs_space, state_space, observations, initial_probs, trans_matrix, emission_matrix)
        trans_matrix, emission_matrix = estimateParametersHMM(observations, hidden_path, state_space_map, obs_space_map)

    return trans_matrix, emission_matrix


def forward_backward(state_space: list[str], obs_space: list[str], observations: str, 
                     initial_probs: npt.NDArray, trans_matrix: npt.NDArray, emission_matrix: npt.NDArray):
    """
    The probability Pr(pi_i = k | x) that the HMM was in state k at step i (for each state k and step i). Pr(pi_i = k, x) is the unconditional probability that a hidden path wil pass through state k at time i and emit x.

    Pr(pi_i = k, x) = sum_{all paths pi with pi_i = k} Pr(x,pi) (This is the forward algorithm k,i)

    sum_{all possible states k, all possible paths} Pr(pi_i=k, x) = 1

    Pr(pi_i = k | x) that the HMM was in state k at step i (for each state k and step i) given it emitted x
    Pr(pi_i = k | x) = Pr(pi_i=k, x) / Pr(x)

    Pr(pi_i = k, x) = (forward algorithm k,i) * backward algorithm k,i

    # https://www.youtube.com/watch?v=yUZ8CBdeJRs
    # https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm
    """
    N = len(state_space)
    K = len(obs_space)
    T = len(observations)

    obs_space_dict = {obs: i for i, obs in enumerate(obs_space)}
    obs_array = [obs_space_dict[obs] for obs in observations]
    obs_array = np.array(obs_array)

    fwd = np.zeros((N, len(observations)), dtype=np.float64)
    alpha_t = initial_probs * emission_matrix[:, obs_array[0]]
    alpha_t = alpha_t / np.linalg.norm(alpha_t, 1)
    fwd[:,0] = alpha_t

    for i, obs in enumerate(obs_array[1:], 1):
        alpha_t = emission_matrix[:, obs] * \
            (np.transpose(trans_matrix) @ alpha_t)
        alpha_t = alpha_t / np.linalg.norm(alpha_t,1)
        fwd[:, i] = alpha_t
    p_fwd = np.sum(alpha_t)

    bkw = np.zeros((N, len(observations)), dtype=np.float64)
    beta_t = np.ones((N,), dtype=np.float64)
    bkw[:, T-1] = beta_t

    for i in range(T-1, 0, -1):
        beta_t = trans_matrix @ np.diag(emission_matrix[:, obs_array[i]]) @ beta_t.transpose()
        beta_t = beta_t / np.linalg.norm(beta_t, 1)
        bkw[:, i-1] = beta_t
    
    p_bkw = np.sum(initial_probs * emission_matrix[:, obs_array[0]] * beta_t)
    a = 0



    print(p_fwd)