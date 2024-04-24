from collections import Counter, defaultdict
from itertools import combinations
import math
import random
import numpy.typing as npt
from typing import Literal
import numpy as np

import networkx as nx
from graphviz import Graph
from BioInfoToolkit.Sequences.Profiling import SequencesProfile, findProfileMostProbableKmers

from BioInfoToolkit.Sequences.StringUtils import getMinimumHammingDistanceToKmer, hamming_distance, kmer_gen, select_random_kmer


def generate_random_sequence(length: int, alphabet: list[str]) -> str:
    seq = "".join(random.choices(alphabet, k=length))
    return seq


def tuple_kmer_gen(seq: tuple, k: int, cyclic: bool = False):
    n = len(seq)

    if k>n or k<0:
        return

    if n == 0:
        yield tuple()
        return

    for i in range(0, n-k+1):
        yield seq[i:i+k]

    if cyclic:
        for i in range(n-k+1, n):
            d = n-i
            yield seq[i:] + seq[:k-d]


def tuple_kmers_frequency_dictionary(sequence: tuple, kmer_len: int, cyclic: bool = False) -> dict[tuple, int]:
    """Builds a dictionary with the frequency of all the k-mers in a sequence

    Args:
        sequence (tuple): tuple sequence
        kmer_len (int): The length of the k-mers

    Returns:
        dict[str, int]: Dictionary with k-mers as keys and the k-mer counts as values
    """

    kmer_frequencies: dict[tuple, int] = dict()

    if kmer_len == 0:
        kmer_frequencies[tuple()] = 1
        return kmer_frequencies

    for kmer in tuple_kmer_gen(sequence, kmer_len, cyclic):
        kmer_frequencies.setdefault(kmer, 0)
        kmer_frequencies[kmer] += 1

    return kmer_frequencies


def complement(seq: str, type: Literal['DNA'] | Literal['RNA']) -> str:
    if type == "DNA":
        mapping = str.maketrans('ATCG', 'TAGC')
    else:
        mapping = str.maketrans('AUCG', 'UAGC')
    return seq.translate(mapping)


def reverseComplement(seq: str, type: Literal['DNA'] | Literal['RNA']) -> str:
    return complement(seq, type)[::-1]


def gc_content(seq: str) -> float:
    """Returns GC constent in a DNA / RNA sequence

    Args:
        seq (str): DNA or RNA string

    Returns:
        _type_: _description_
    """
    return (seq.count('C') + seq.count('G')) / (len(seq))


def minimum_skew(sequence: str):
    """Define the skew of a DNA string Genome, denoted Skew(Genome), as the difference between the total number of occurrences of 'G' and 'C' in Genome. Let Prefixi (Genome) denote the prefix (i.e., initial substring) of Genome of length i.

    Args:
        sequence (str): _description_

    Returns:
        _type_: _description_
    """
    skew_val = 0
    skews: list[int] = [skew_val]
    for char in sequence:
        if char == 'G':
            skew_val += 1
        elif char == 'C':
            skew_val -= 1
        skews.append(skew_val)

    min_skew = min(skews)
    positions = [i for i, skew_val in enumerate(skews) if skew_val == min_skew]
    return min_skew, positions


def distanceMatrix(seqs: list[str]) -> npt.NDArray[np.float64]:
    n = len(seqs)
    distMat = np.zeros((n, n), dtype=np.float64)
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs[i + 1:], i + 1):
            dist = hamming_distance(seq1, seq2) / min(len(seq1), len(seq2))
            distMat[i, j] = dist
            distMat[j, i] = dist
    return distMat


def count_transitions_transversions(seq1: str, seq2: str) -> tuple[int, int]:
    """ Counts transitions and transversions in 2 sequences of the same length
        Transitions A <-> G, C<->T
        Transversions A<->C, A<->T, G<->C, G<->T

    Args:
        seq1 (str): Sequence 1
        seq2 (str): Sequence 2

    Returns:
        tuple[int, int]: (transitions, transversions)
    """
    AG = ("A", "G")
    CT = ("C", "T")
    transitions = 0
    transversions = 0
    for a, b in zip(seq1, seq2):
        if a == b:
            continue
        if (a in AG and b in AG) or (a in CT and b in CT):
            transitions += 1
        else:
            transversions += 1

    return transitions, transversions


def failure_array(seq: str) -> list[int]:
    """Returns the failure array of a sequence (seq):
        The failure array of s is an array P of length n for which P[k] is the length of the longest substring
        s[j:k+1] that is equal to some prefix s[:k−j+1], where j cannot equal 0
        (otherwise, P[k] would always equal k). By convention, P[0]=0.

    Args:
        seq (str): Sequence

    Returns:
        list[int]: Failure array
    """
    n = len(seq)
    failure: list[int] = []

    if n == 0:
        return failure
    failure.append(0)
    for k in range(1, n):
        pk = 0
        pk_1 = failure[k - 1]  # prev Pk
        for j in range(k - pk_1, k + 1):
            prefix = seq[:k - j + 1]
            substring = seq[j:k + 1]
            if prefix == substring:
                pk = len(prefix)
                break
        failure.append(pk)
    return failure


def is_subsequence_str(seq: str, subseq: str):
    i1 = 0
    v1 = subseq[i1]
    l = len(subseq)
    for v2 in seq:
        if v2 != v1:
            continue
        i1 += 1
        if i1 == l:
            return True
        v1 = subseq[i1]
    return False


def longestCommonString(sequences: list[str]) -> str:
    """Find the longest common substring from a list of strings

    Args:
        sequences (list[str]): string list

    Returns:
        str: longest common substring
    """
    n = len(sequences)
    seq0 = sequences[0]
    l = len(seq0)
    res = ''

    for i in range(l):
        for j in range(i+1, l):
            substr = seq0[i:j]
            is_common_to_all = all(substr in seq for seq in sequences[1:])
            if is_common_to_all and len(substr) > len(res):
                res = substr

    return res


def longestCommonSubsequenceTable(seq1: str, seq2: str) -> list[list[int]]:
    """Computes the longest common subsequence matrix

    https://en.wikipedia.org/wiki/Longest_common_subsequence

    Args:
        seq1 (str):
        seq2 (str):

    Returns:
        list[list[int]]: longest common subsequence matrix C[i][j] contains the length 
        of longest common subsequence for seq1[0:i] and seq2[0:j]
    """

    m = len(seq1)
    n = len(seq2)

    C: list[list[int]] = [[0 for _ in range(n+1)] for _ in range(m+1)]

    for i, xi in enumerate(seq1, 1):
        for j, yj in enumerate(seq2, 1):
            if xi == yj:
                C[i][j] = C[i-1][j-1] + 1
            else:
                C[i][j] = max(C[i][j-1], C[i-1][j])

    return C


def longestCommonSubsequenceBacktrack(C: list[list[int]],
                                      seq1: str, seq2: str,
                                      i: int, j: int) -> str:
    """Computes the longest common subsequence (LCS), iteratively.
    To compute the LCS between seq1 and seq2 call this function with i=len(seq1), j=(seq2),
    C = computeLongestCommonSubsequenceTable(seq1, seq2)

    Args:
        C (list[list[int]]): longest common subsequence matrix
        seq1 (str):
        seq2 (str):
        i (int): seq1 index to calculate the LCS from
        j (int): seq2 index to calculate the LCS from

    Returns:
        str: longest common substring
    """
    result = ""
    while i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]:
            result = seq1[i-1] + result
            i, j = i - 1, j - 1
        elif C[i][j-1] > C[i-1][j]:
            j = j - 1
        else:
            i = i - 1
    return result


def longestCommonSubsequence(seq1: str, seq2: str):
    """Computes the longest common suqsequence between two strings
    https://en.wikipedia.org/wiki/Longest_common_subsequence

    Args:
        seq1 (str): 
        seq2 (str): 

    Returns:
        str: longest common substring
    """

    C = longestCommonSubsequenceTable(seq1, seq2)
    lcs = longestCommonSubsequenceBacktrack(
        C, seq1, seq2, len(seq1), len(seq2))
    return lcs, C


def shortestCommonSupersequence(s: str, t: str) -> str:
    m, n = len(s), len(t)

    # Initialize a table to store the lengths of the shortest common supersequences
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Fill in the table using dynamic programming
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j
            elif j == 0:
                dp[i][j] = i
            elif s[i - 1] == t[j - 1]:
                dp[i][j] = 1 + dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i - 1][j], dp[i][j - 1])

    # Reconstruct the shortest common supersequence
    result = []
    i, j = m, n
    while i > 0 and j > 0:
        if s[i - 1] == t[j - 1]:
            result.append(s[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] < dp[i][j - 1]:
            result.append(s[i - 1])
            i -= 1
        else:
            result.append(t[j - 1])
            j -= 1

    # Add remaining characters from s and t
    while i > 0:
        result.append(s[i - 1])
        i -= 1
    while j > 0:
        result.append(t[j - 1])
        j -= 1

    # Reverse the result to get the correct order
    result.reverse()

    return ''.join(result)


def LongestCommonSubsequenceBacktrackAll(C: list[list[int]],
                                         seq1: str, seq2: str,
                                         i: int, j: int) -> set[str]:
    if i == 0 or j == 0:
        return {""}
    if seq1[i-1] == seq2[j-1]:
        return {x + seq1[i-1] for x in LongestCommonSubsequenceBacktrackAll(C, seq1, seq2, i-1, j-1)}

    R: set[str] = set()
    if C[i][j-1] >= C[i-1][j]:
        R = LongestCommonSubsequenceBacktrackAll(C, seq1, seq2, i, j-1)
    if C[i-1][j] >= C[i][j-1]:
        R = R.union(LongestCommonSubsequenceBacktrackAll(
            C, seq1, seq2, i-1, j))
    return R


def bondingGraph(seq: str, bonding_map: dict[str, list[str]],
                 min_dist: int = 1, add_adjacencies: bool = False) -> nx.MultiGraph:
    """ Creates a bonding graph for a given sequence.
    Given a string s=s1…sn, a bonding graph for s is formed as follows. First, assign each symbol of s to a node, and arrange these nodes in order around a circle, connecting them with edges called adjacency edges. Second, form all possible edges allowed by the bonding map, called basepair edges; we will represent basepair edges with dashed edges, as illustrated by the bonding graph. 

    For an RNA string the possible edges are {A, U} and {C, G}.

    Args:
        seq (str): _description_
        bonding_map (dict[str, list[str]], optional): _description_.
        min_dist (int, optional): _description_. Defaults to 1.
        add_adjacencies (bool, optional): Wether to add adjacency edges to the graph. Defaults to False.

    Returns:
        nx.MultiGraph: Note that if adjacency edges are added there might be 2 edges between some pair of adjacent nodes (a MultiGraph is necessary)
    """
    n = len(seq)
    gr = nx.MultiGraph()
    gr.add_nodes_from([(i, {"nt": nt}) for i, nt in enumerate(seq)])

    # adjacency edges
    if add_adjacencies:
        gr.add_edges_from(
            [(i, i + 1, {"type": "adjacency"}) for i in range(n - 1)])

    # bonding edges
    for i, j in combinations(range(n), 2):
        diff = abs(i-j)
        dist = min(diff, n-diff)
        if seq[j] in bonding_map[seq[i]] and dist >= min_dist:
            gr.add_edge(i, j, type="bonding")

    return gr


def drawBondingGraph(g: nx.MultiGraph) -> Graph:
    """Drawns a bonding graph using the neato engine

    Args:
        g (nx.MultiGraph): the bonding graph (networkx MultiGraph)

    Returns:
        Graph: drawn graph
    """
    n = g.number_of_nodes()
    r = 0.4*n

    def getx(i: int) -> float:
        return r * math.cos(2 * math.pi * i / n)

    def gety(i: int) -> float:
        return r * math.sin(2 * math.pi * i / n)

    node_pos = {node: ",".join(
        [str(getx(node)), str(gety(node))])+'!' for node in g.nodes}

    neato = Graph(format='svg', engine="neato")
    for node, nattr in g.nodes.items():
        nt = nattr['nt']
        neato.node(str(node), f"{nt}", pos=node_pos[node])
    for u, v, d in g.edges(data=True):
        if d['type'] == "adjacency":
            neato.edge(str(u), str(v), style="bold")
        else:
            neato.edge(str(u), str(v), style="dashed")
    return neato


def countNoncrossingMatches(seq: str,
                            bonding_map: dict[str, list[str]],
                            min_dist: int = 1,
                            mod: int | None = None,
                            cache: dict[str, int] = dict()) -> int:
    """Counts the total number of (not necessarily perfect) noncrossing matches in a given sequence

    Args:
        seq (str): The sequence to be analized
        bonding_map (dict[str, list[str]]): Nucleotide bonding map
        min_dist (int): The minimum distance between nucleotides (mod len(seq)) 
            necessary for the nucleotides to bond
        cache (dict[str, int], optional): Memoization cache. Saves the values of computed sequences. 
            Defaults to dict().

    Returns:
        int: Total count of noncrossing matches
    """
    if seq in cache:
        return cache[seq]

    motzkinNum = 0
    n = len(seq)
    if n == 0:
        motzkinNum = 1
        cache[seq] = motzkinNum
        return motzkinNum

    node1 = 0
    nt1 = seq[0]
    adjNodes = [i for i, nt2 in enumerate(
        seq[1:], 1) if nt2 in bonding_map[nt1]]

    seq0 = seq[1:]
    c0 = countNoncrossingMatches(seq0, bonding_map, min_dist, mod, cache)

    motzkinNum = (motzkinNum + c0) if mod is None else (motzkinNum + c0) % mod

    for node2 in adjNodes:
        if abs(node1-node2) < min_dist:
            continue
        seq1 = seq[node1+1:node2]
        seq2 = seq[node2+1:]

        c1 = countNoncrossingMatches(
            seq1, bonding_map, min_dist, mod, cache)
        c2 = countNoncrossingMatches(
            seq2, bonding_map, min_dist, mod, cache)

        aux = (motzkinNum + c1 * c2)
        motzkinNum = aux if mod is None else aux % mod

    cache[seq] = motzkinNum

    return motzkinNum


def countNoncrossingPerfectMatches(seq: str, bonding_map: dict[str, list[str]], 
                                   min_dist: int = 1, mod: int = 1000000, 
                                   cache: dict[str, int] = dict()) -> int:
    if seq in cache:
        return cache[seq]
    
    def canBeMatchedByCount(seq: str) -> bool:
        counts = Counter(seq)
        return counts.get('A', 0) == counts.get('U', 0) and counts.get('G', 0) == counts.get('C', 0)

    catalanNum = 0
    n = len(seq) 
    if n == 0:
        catalanNum = 1
        cache[seq] = catalanNum
        return catalanNum
    if n % 2 != 0:
        cache[seq] = catalanNum
        return catalanNum
    if not canBeMatchedByCount(seq):
        cache[seq] = catalanNum
        return catalanNum

    node1 = 0
    nt1 = seq[0]
    # nodes that can be connected to the first nucleotide/symbol allowed by the bonding map
    adjNodes = [i for i, nt2 in enumerate(seq[1:], 1) if nt2 in bonding_map[nt1]]

    for node2 in adjNodes:
        if abs(node1-node2) < min_dist:
            continue
        seq1 = seq[node1+1:node2]
        seq2 = seq[node2+1:]

        c1 = countNoncrossingPerfectMatches(seq1, bonding_map, min_dist, mod, cache)
        c2 = countNoncrossingPerfectMatches(seq2, bonding_map, min_dist, mod, cache)

        catalanNum = (catalanNum + c1 * c2) % mod
    cache[seq] = catalanNum

    return catalanNum


def deBruijnGraph(seqs: list[str], k: int) -> defaultdict[str, set[str]]:
    S = set(seqs)
    Src = set(reverseComplement(seq, 'DNA') for seq in S)
    SSrc = S.union(Src)

    graph: defaultdict[str, set[str]] = defaultdict(set)
    for k1mer in SSrc:
        l = len(k1mer)
        n1 = k1mer[:k]
        n2 = k1mer[l-k:]
        graph[n1].add(n2)
    return graph


def greedyMotifSearch(sequences: list[str], k: int, alphabet: list[str], pseudocounts: int = 0):
    """Implements the greedy motif search algorithm with pseudocounts

    Args:
        sequences (list[str]): sequence to search
        k (int): motif length
        alphabet (list[str]): _description_
        pseudocounts (int, optional): Defaults to 0 (no pseudocounts).

    Returns:
        best_motifs, best_score (tuple[list[str], int]): 
    """
    bestMotifs = [seq[0:k] for seq in sequences]
    scoreBestMotifs = sum(sum(getMinimumHammingDistanceToKmer(seq, kmer)
                          for seq in sequences) for kmer in bestMotifs)

    for kmer in kmer_gen(sequences[0], k):
        motifs = [kmer]

        for _, seq2 in enumerate(sequences[1:]):
            profile = SequencesProfile(motifs, alphabet, pseudocounts)
            prob_mat: list[list[float]
                           ] = profile.get_probability_matrix().tolist()
            kmers, _ = findProfileMostProbableKmers(
                seq2, k, prob_mat, alphabet)
            motif = kmers[0]
            motifs.append(motif)

        score = sum(sum(getMinimumHammingDistanceToKmer(seq, kmer)
                        for seq in sequences) for kmer in motifs)
        if score < scoreBestMotifs:
            scoreBestMotifs = score
            bestMotifs = motifs
    return bestMotifs, scoreBestMotifs


def randomizedMotifSearchMonteCarlo(sequences: list[str], k: int, alphabet: list[str], n_iter: int = 1000, pseudocounts: int = 0):
    bestMotifs: list[str] = []
    bestScore = math.inf

    for i in range(n_iter):
        motifs, score = randomizedMotifSearch(
            sequences, k, alphabet, pseudocounts)
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
    return bestMotifs, bestScore


def randomizedMotifSearch(sequences: list[str], k: int, alphabet: list[str], pseudocounts: int = 0):
    motifs = [select_random_kmer(seq, k) for seq in sequences]
    bestMotifs = motifs
    bestProfile = SequencesProfile(bestMotifs, alphabet, pseudocounts)
    scoreBestMotifs = bestProfile.score()

    while True:
        profile = SequencesProfile(motifs, alphabet, pseudocounts)

        prob_mat: list[list[float]] = profile.get_probability_matrix().tolist()
        motifs: list[str] = []
        for seq in sequences:
            kmers, _ = findProfileMostProbableKmers(
                seq, k, prob_mat, alphabet)
            motif = kmers[0]
            motifs.append(motif)

        profile = SequencesProfile(motifs, alphabet, pseudocounts)
        score = profile.score()
        if score < scoreBestMotifs:
            scoreBestMotifs = score
            bestMotifs = motifs
        else:
            return bestMotifs, scoreBestMotifs


def gibbsSamplerMonteCarlo(sequences: list[str], k: int, alphabet: list[str], N: int,
                           n_starts: int, pseudocounts: int = 1) -> tuple[list[str], int|float]:
    """Runs Monte Carlo search using the gibbs sampler

    Args:
        sequences (list[str]): sequences to search
        k (int): motif length
        alphabet (list[str]):
        N (int): number of iterations of the sampler
        n_starts (int, optional): Number of random starts (runs the sampler N times).
        pseudocounts (int, optional): Defaults to 1.

    Returns:
        best_motifs, best_score (tuple[list[str], int]):
    """
    best_motifs: list[str] = []
    best_score = math.inf

    for _ in range(n_starts):
        motifs, score = gibbsSampler(sequences, k, N, alphabet, pseudocounts)
        if score < best_score:
            best_score = score
            best_motifs = motifs

    return best_motifs, best_score


def gibbsSampler(sequences: list[str], k: int, N: int, alphabet: list[str], pseudocounts: int = 1):
    motifs = [select_random_kmer(seq, k) for seq in sequences]
    bestMotifs = motifs
    profile = SequencesProfile(bestMotifs, alphabet, pseudocounts)
    scoreBestMotifs = profile.score()
    t = len(sequences)

    alphabet_map = {symb: i for i, symb in enumerate(alphabet)}

    for _ in range(N):
        i = random.randrange(0, t)
        motifs2 = [motif for i2, motif in enumerate(motifs) if i2 != i]
        profile = SequencesProfile(motifs2, alphabet, pseudocounts)
        prob_mat: list[list[float]] = profile.get_probability_matrix().tolist()

        probs = [math.prod(prob_mat[alphabet_map[symb]][i]
                           for i, symb in enumerate(kmer)) for kmer in kmer_gen(sequences[i], k)]

        idx = random.choices(
            list(range(len(sequences[i])-k+1)), k=1, weights=probs)[0]
        kmer = sequences[i][idx:idx+k]
        motif_i = kmer
        motifs[i] = motif_i

        profile = SequencesProfile(motifs, alphabet, pseudocounts)
        score = profile.score()
        if score < scoreBestMotifs:
            scoreBestMotifs = score
            bestMotifs = motifs

    return bestMotifs, scoreBestMotifs


def distanceBetweenPatternAndStrings(strings: list[str], pattern: str):
    k = len(pattern)

    distance = 0
    for string in strings:
        dist2 = math.inf
        for kmer in kmer_gen(string, k):
            dist3 = hamming_distance(pattern, kmer)
            if dist3 < dist2:
                dist2 = dist3
        distance += dist2
    return distance

