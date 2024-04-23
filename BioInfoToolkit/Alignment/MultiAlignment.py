

from itertools import combinations, product
from typing import Callable
import numpy as np


def getMultiAlignScoringFunc(matchScore: int = 1, mismatchScore: int = -1,
                             similarityDict: dict[tuple[str, str], int] | None = None) -> Callable[[list[str]], int]:
    def score(chars: list[str]) -> int:
        totalScore = 0

        for x, y in combinations(chars, 2):
            if similarityDict:
                totalScore = similarityDict[(x, y)]
            else:
                similarityScore = matchScore if x == y else mismatchScore
                totalScore += similarityScore

        return totalScore

    return score


def multiAlignmentCost(aligned_seqs: list[str], scoringFunc: Callable[[list[str]], int]):
    score = 0
    for symbs in zip(*aligned_seqs):
        score += scoringFunc(list(symbs))
    return score


def multisequenceAlign(seqs: list[str],
                       scoringFunc: Callable[[list[str]], int]) -> tuple[list[str], int]:
    """Needleman-Wunsch algorithm to compute optimal global sequence alignment. Extends the 2 sequence alignement algorithm to arbitrary n sequences. Creates a matrix of size (n1, n2, ..., nk), where k is the number of sequences and ni is the length of sequence i

    Args:
        seqs (list[str]): a list of sequences to align

    Returns:
        tuple[list[str], int]: tuple with a list of aligned sequences and the alignment score
    """

    # num_seq = len(seqs)
    seqs_len = tuple(len(seq) for seq in seqs)

    def tuple_add(t1: tuple[int, ...], t2: tuple[int, ...]):
        return tuple(a + b for a, b in zip(t1, t2))

    def generate_predecessors_delta(idxs: tuple[int, ...]):
        aux = tuple([-1, 0] if i > 0 else [0] for i in idxs)
        for _, di in enumerate(product(*aux)):
            if all(v == 0 for v in di):
                continue
            yield di

    # score matrix
    fDims = tuple(s+1 for s in seqs_len)
    F = np.zeros(fDims, dtype=np.int32)

    for idxs in product(*[list(range(0, n+1)) for n in seqs_len]):
        values: list[int] = []
        for delta in generate_predecessors_delta(idxs):
            chars = ['-' if k == 0 else seqs[j][idxs[j] + k]
                     for j, k in enumerate(delta)]
            value = F[tuple_add(idxs, delta)]+scoringFunc(chars)
            values.append(value)

        max_val = max(values) if len(values) else 0
        F[idxs] = max_val

    # backtrack
    alignedSeqs = ['' for _ in seqs]
    idxs = seqs_len
    final_score = int(F[idxs])

    while any(i > 0 for i in idxs):
        # find which delta we came from
        for delta in generate_predecessors_delta(idxs):
            chars = ['-' if k == 0 else seqs[j][idxs[j] + k]
                     for j, k in enumerate(delta)]
            if F[idxs] == F[tuple_add(idxs, delta)]+scoringFunc(chars):
                alignedSeqs = [chars[j] + seq for j,
                               seq in enumerate(alignedSeqs)]
                idxs = tuple_add(idxs, delta)
                break

    return alignedSeqs, final_score
