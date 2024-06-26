from collections import Counter, OrderedDict
import math
import numpy as np
import numpy.typing as npt

from BioInfoToolkit.Sequences.StringUtils import kmer_gen


def build_profile_count_matrix(sequences: list[str], alphabet: list[str], pseudocounts: int = 0):
    n = len(sequences[0])
    k = len(alphabet)
    count_matrix = np.full((k,n), pseudocounts, dtype=np.uint64)

    alphabet_map = {c: i for i, c in enumerate(alphabet)}
    for seq in sequences:
        for i, c in enumerate(seq):
            j = alphabet_map[c]
            count_matrix[j][i] += 1
    return count_matrix


class SequencesProfile:
    """Creates a profile of the sequences strings. The length of the profile is equal to the minimum length of a
        string in sequences. Example:

        DNA Strings
            A T C C A G C T
            G G G C A A C T
            A T G G A T C T
            A A G C A A C C
            T T G G A A C T
            A T G C C A T T
            A T G G C A C T

        Profile
            A   5 1 0 0 5 5 0 0
            C   0 0 1 4 2 0 6 1
            G   1 1 6 3 0 1 0 0
            T   1 5 0 0 0 1 1 6
    """
    alphabet: list[str]
    count_matrix: npt.NDArray[np.uint64]

    def __init__(self, sequences: list[str], alphabet: list[str] | None, pseudocounts: int = 0) -> None:
        used_symbols = sorted(list(set(symbol for seq in sequences for symbol in seq)))
        if not alphabet:
            alphabet = used_symbols
        else:
            alphabet = sorted(set(used_symbols).union(set(alphabet)))
        self.alphabet = alphabet

        count_matrix = build_profile_count_matrix(sequences, alphabet, pseudocounts)
        self.count_matrix = count_matrix
    
    def __repr__(self) -> str:
        out = ""
        n = len(self.count_matrix[0])
        lengths = [max((len(str(v[i])) for v in self.count_matrix)) for i in range(n)]
        for i, symbol in enumerate(self.alphabet):
            values = [str(v).rjust(lengths[j]) for j, v in enumerate(self.count_matrix[i])]
            line = f"{symbol}: " + " ".join(values) + "\n"
            out += line
        return out.removesuffix('\n')
    
    def __str__(self) -> str:
        return self.__repr__()

    @classmethod
    def from_profile_matrix(cls, alphabet: list[str], profile_mat: list[list[int]]) -> "SequencesProfile":
        n = len(profile_mat)
        seqs = ['' for _ in range(n)]
        profile = SequencesProfile(seqs, alphabet)
        count_matrix = np.array(profile_mat, dtype=np.uint64)
        profile.count_matrix = count_matrix
        return profile

    def consensus(self) -> str:
        """Returns a possible consensus, if multiple consensus are possible it will only return one

        Returns:
            str: consensus string
        """
        seq = ""
        max_counts = np.max(self.count_matrix, 0)
        for i, max_val in enumerate(max_counts):
            nucs = [self.alphabet[j] for j,v in enumerate(self.count_matrix[:,i]) if v == max_val]
            seq += nucs[0]
        return seq
    
    def possible_consensus(self) -> list[list[str]]:
        consensus: list[list[str]] = []
        max_counts = np.max(self.count_matrix, 0)
        for i, max_val in enumerate(max_counts):
            nucs = [self.alphabet[j] for j, v in enumerate(
                self.count_matrix[:, i]) if v == max_val]
            consensus.append(nucs)
        return consensus
    
    def is_consensus(self, consensus: str) -> bool:
        max_counts = np.max(self.count_matrix, 0)
        alphabet_map = {c: i for i, c in enumerate(self.alphabet)}
        for i, max_val in enumerate(max_counts):
            j = alphabet_map[consensus[i]]
            if self.count_matrix[j,i] != max_val:
                return False
        return True

    def get_probability_matrix(self) -> npt.NDArray[np.float64]:
        prob_matrix = np.array(self.count_matrix, dtype=np.float64)
        sums = np.sum(prob_matrix, 0)

        prob_matrix = prob_matrix / sums
        return prob_matrix


    def to_list_matrix(self) -> list[list[int]]:
        profile_mat: list[list[int]] = self.count_matrix.tolist()
        return profile_mat


    def score(self) -> int:
        sums = np.sum(self.count_matrix, 0)
        max_counts = np.max(self.count_matrix, 0)
        score = int(np.sum(sums-max_counts))
        return score


def findProfileMostProbableKmers(text: str, k: int, profile_mat: list[list[float]], 
                                 alphabet: list[str]) -> tuple[list[str], float]:
    """Given a profile matrix Profile, we can evaluate the probability of every k-mer in a string Text and find a Profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated by Profile among all k-mers in Text.

    Args:
        text (str): a string
        k (int): k-mer length
        profile_mat (list[list[float]]): an n x k matrix of probabilities
        alphabet (set[str]): An alphabet of length n

    Returns:
        most_prob_strings, prob (tuple[list[str], float]): tuple with a list of the most probable strings and corresponding probability
    """
    alphabet = sorted(alphabet)
    alphabet_map = {symb: i for i, symb in enumerate(alphabet)}

    best_prob = float(0)
    kmers: list[str] = []

    for kmer in kmer_gen(text, k):
        prob = math.prod(profile_mat[alphabet_map[symb]][i]
                         for i, symb in enumerate(kmer))
        if prob > best_prob:
            best_prob = prob
            kmers = [kmer]
        elif prob == best_prob and kmer not in kmers:
            kmers.append(kmer)

    return kmers, best_prob
