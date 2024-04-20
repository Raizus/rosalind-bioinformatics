from collections import defaultdict
from BioInfoToolkit.IO import readTextFile
from BioInfoToolkit.Alignment.Alignment import SimilarityScore, add_backtrack_paths_to_alignment_svg_global_alignement, drawAlignmentMatrixSvg, fittingAlignmentLinearGap
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import time 

"""
https://rosalind.info/problems/ksim/

    Given: A positive integer k (k≤50), a DNA string s of length at most 5 kbp representing a motif, and a DNA string t of length at most 50 kbp representing a genome.

    Return: All substrings t' of t such that the edit distance dE(s,t') is less than or equal to k. Each substring should be encoded by a pair containing its location in t followed by its length.

"""

OutputT = list[tuple[int, int]]


def findSimilarMotifs(seq: str, motif: str, max_dist: int):
    """Finds all similiar motifs to motif that align with a substring of seq, i.e. it finds all substrings seq' of seq such that the edit distance edit_distance(seq', motif) <= max_dist

    Args:
        seq (str): _description_
        motif (str): _description_
        de (int): _description_
    """

    gap_penalty = 1
    matchScore = 0
    mismatchScore = 1
    similarity_score = SimilarityScore(matchScore, mismatchScore)

    m = len(seq)
    n = len(motif)

    # observations:
    # since we're working with long sequences, and the cost of alignment always increases, we can reduce our search space by pruning coordinates whose cost as increased above max_dist
    # for every column we must keep track of all costs as long as they're below max_dist, and origin coordinates that
    # lead to it

    # defaultdict[int, set[int]] contains the alignment cost and the corresponding idx of seq where the alignment(s) start
    # the key's of the dictionary correspond to the index of seq and costs associated for the last column
    last_col_dict: dict[int, defaultdict[int, set[int]]] = dict()
    cost_dict = defaultdict(set)
    cost_dict[0] = set([0])
    last_col_dict[0] = cost_dict
    for i in range(1, m+1):
        cost_dict = defaultdict(set)
        cost_dict[0] = set([i])
        last_col_dict[i] = cost_dict
        up = last_col_dict.get(i-1, defaultdict(set))
        for cost, from_i in up.items():
            insert = cost + gap_penalty
            if insert <= max_dist:
                cost_dir = last_col_dict.setdefault(i, defaultdict(set))
                cost_dir[insert].update(from_i)

    for j, tj in enumerate(motif, 1):
        curr_col_dict: dict[int, defaultdict[int, set[int]]]  = dict()
        if j*gap_penalty <= max_dist:
            curr_col_dict[0] = defaultdict(set)
            curr_col_dict[0][j*gap_penalty] = set([0])

        for i, si in enumerate(seq, 1):
            diag = last_col_dict.get(i-1, defaultdict(set))
            for cost, from_i in diag.items():
                score = similarity_score.score(si, tj)
                match = cost + score
                if match <= max_dist:
                    cost_dir = curr_col_dict.setdefault(i, defaultdict(set))
                    cost_dir[match].update(from_i)
                
            left = last_col_dict.get(i, defaultdict(set))
            for cost, from_i in left.items():
                delete = cost + gap_penalty
                if delete <= max_dist:
                    cost_dir = curr_col_dict.setdefault(i, defaultdict(set))
                    cost_dir[delete].update(from_i)
            
            up = curr_col_dict.get(i-1, defaultdict(set))
            for cost, from_i in up.items():
                insert = cost + gap_penalty
                if insert <= max_dist:
                    cost_dir = curr_col_dict.setdefault(i, defaultdict(set))
                    cost_dir[insert].update(from_i)

        last_col_dict = curr_col_dict

    idx_pairs: OutputT = []
    for end, cost_dict in curr_col_dict.items():
        for cost, from_i in cost_dict.items():
            for start in from_i:
                pair = (start, end)
                idx_pairs.append(pair)

    return idx_pairs


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(seq1: str, seq2: str, k: int) -> OutputT:
    idx_pairs = findSimilarMotifs(seq1, seq2, k)
    start_length_array: list[tuple[int,int]] = []

    for start, end in idx_pairs:
        start_length_array.append((start+1, end-start))
    return start_length_array


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    start_length_array: OutputT = []
    for line in lines:
        i1, i2 = [int(v) for v in line.split()]
        start_length_array.append((i1, i2))

    return start_length_array


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    k = int(lines[0])
    s = lines[1]
    t = lines[2]

    result = solve(t, s, k)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ksim_2.txt'

    lines = readTextFile(path)
    k = int(lines[0])
    s = lines[1]
    t = lines[2]

    # print(len(t))
    # print(len(s))
    # print()

    # t1 = time.time()
    start_length_array = solve(t, s, k)
    # t2 = time.time()
    # print(t2-t1)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for pair in start_length_array:
        out = ' '.join(str(p) for p in pair)
        print(out)
        writeTextFile(result_path, out, 'a')

    # correct = solve_and_check(path)
    # print(correct)

    # gap_penalty = -1
    # matchScore = 0
    # mismatchScore = -1
    # similarity_score = SimilarityScore(matchScore, mismatchScore)
    # H, originMat, _ = fittingAlignmentLinearGap(
    #     t, s, gap_penalty, similarity_score)

    # idxs = [(i, len(s)) for i in range(len(t)+1) if H[i][len(s)] >= -k]
    # d = drawAlignmentMatrixSvg(H, t, s)
    # add_backtrack_paths_to_alignment_svg_global_alignement(
    #     d, originMat, t, s, idxs)
    # d.save_svg('ksim.svg')
