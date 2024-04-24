from collections import defaultdict
from itertools import combinations, product
import math
import random
import re


def strip_inserts(seq: str) -> str:
    """Replaces all '-' in a string with emptry strings ''

    Args:
        seq (str):

    Returns:
        str:
    """
    return seq.replace('-', '')


def rotationally_equivalent(s1: str, s2: str) -> bool:
    res = False
    for idx in range(len(s1)):
        if s1[idx:] + s1[:idx] == s2:
            res = True
            break
    return res


def patternToNumber(string: str, alphabet: list[str]):
    L = len(alphabet)
    convert_dict = {symb: i for i, symb in enumerate(alphabet)}

    def _patternToNumber(string: str):
        if len(string) == 0:
            return 0
        symbol = string[-1]
        return L * _patternToNumber(string[:-1]) + convert_dict[symbol]
    
    return _patternToNumber(string)


def numberToPattern(index: int, k: int, alphabet: list[str]):
    L = len(alphabet)
    convert_dict = {i: symb for i, symb in enumerate(alphabet)}

    def _numberToPattern(ind: int, k: int):
        if k == 1:
            return convert_dict[ind]
        prefix_index = ind // L
        r = ind % L
        symbol = convert_dict[r]
        prefix_pattern = _numberToPattern(prefix_index, k-1)
        return prefix_pattern + symbol
    
    return _numberToPattern(index, k)


def select_random_kmer(seq: str, k: int) -> str:
    idx = random.randrange(0, len(seq)-k+1)
    return seq[idx:idx+k]


def kmer_gen(seq: str, k: int, cyclic: bool = False):
    """Generates k-mers of length k from a sequence

    Args:
        seq (str): _description_
        k (int):
        cyclic (bool, optional): Indicates if the sequence is cyclic. Defaults to False.

    Yields:
        _type_: _description_
    """
    n = len(seq)

    if k > n or k < 0:
        return

    if n == 0:
        yield ''
        return

    for i in range(0, n-k+1):
        yield seq[i:i+k]

    if cyclic:
        for i in range(n-k+1, n):
            d = n-i
            yield seq[i:] + seq[:k-d]


def kmer_composition(seq: str, k: int, cyclic: bool = False) -> list[str]:
    """Returns the kmer composition of a sequence

    Args:
        seq (str): 
        k (int): 
        cyclic (bool, optional): _description_. Defaults to False.

    Returns:
        set[str]: kmer composition
    """
    kmers: set[str] = set()
    for kmer in kmer_gen(seq, k, cyclic):
        kmers.add(kmer)
    return sorted(kmers)


def generate_kmers_from_alphabet(k: int, alphabet: list[str]):
    """Generates all possible k-mers of length k given an alphabet

    Args:
        k (int): _description_
        alphabet (set[str]): _description_

    Yields:
        str: kmer
    """
    alphabet = sorted(alphabet)
    for symbs in product(alphabet, repeat=k):
        kmer = ''.join(symbs)
        yield kmer


def get_kmers(seq: str, k: int, cyclic: bool = False) -> list[str]:
    """Given a string and an integer k, returns a list of all k-mers in that string

    Args:
        seq (str): A sequence
        k (int): k-mer length

    Returns:
        list[str]: list of k-mers
    """
    kmers: list[str] = [kmer for kmer in kmer_gen(seq, k, cyclic)]
    return kmers


def count_kmer(sequence: str, kmer: str) -> int:
    """Counts repeating k-mers in a sequence, including overlapping k-mers

    Args:
        sequence (str):
        kmer (str): target k-mer

    Returns:
        int: count
    """

    k = len(kmer)
    kmer_count = 0
    for kmer2 in kmer_gen(sequence, k):
        if kmer2 == kmer:
            kmer_count += 1

    return kmer_count


def kmers_frequency_dictionary(sequence: str, kmer_len: int, cyclic: bool = False) -> dict[str, int]:
    """Builds a dictionary with the frequency of all the k-mers in a sequence

    Args:
        sequence (str): The bio sequence
        kmer_len (int): The length of the k-mers

    Returns:
        dict[str, int]: Dictionary with k-mers as keys and the k-mer counts as values
    """

    kmer_frequencies: dict[str, int] = dict()

    if kmer_len == 0:
        kmer_frequencies[''] = 1
        return kmer_frequencies

    for kmer in kmer_gen(sequence, kmer_len, cyclic):
        kmer_frequencies.setdefault(kmer, 0)
        kmer_frequencies[kmer] += 1

    return kmer_frequencies


def kmers_frequency_array(text: str, k: int, alphabet: list[str]) -> list[int]:
    kmers_dict = kmers_frequency_dictionary(text, k)
    array: list[int] = []
    for kmer in generate_kmers_from_alphabet(k, alphabet):
        array.append(kmers_dict.get(kmer, 0))
    return array


def find_most_frequent_kmers(string: str, kmer_len: int) -> list[str]:
    """Finds the most frequent k-mers of length kmer_len

    Args:
        sequence (str): The string to search
        kmer_len (int): The length of the k-mers to search for

    Returns:
        list[str]: List of the most frequent k-mers
    """
    kmer_frequencies: dict[str, int] = kmers_frequency_dictionary(
        string, kmer_len)

    highest_freq = max(kmer_frequencies.values())
    kmers = [kmer for kmer, freq in kmer_frequencies.items() if freq ==
             highest_freq]
    return kmers


def paired_kmers_gen(sequence: str, k: int, d: int):
    """Generates the paired composition of a sequence

    Args:
        sequence (str): _description_
        k (int): k-mer length
        d (int): distance between the start of each pair

    Yields:
        _type_: _description_
    """
    n = len(sequence)
    m = n - (d + k) + 1
    for i in range(m):
        kmer1 = sequence[i:i+k]
        kmer2 = sequence[i+d:i+d+k]
        yield (kmer1, kmer2)


def paired_composition(sequence: str, k: int, d: int) -> list[tuple[str, str]]:
    """Returns the paired composition of a sequence, in lexicographic order

    Args:
        sequence (str):
        k (int): k-mer length
        d (int): distance between the start of each pair

    Returns:
        list[tuple[str, str]]: paired composition, sorted lexicographically
    """
    comp: set[tuple[str, str]] = set()
    for paired_kmer in paired_kmers_gen(sequence, k, d):
        comp.add(paired_kmer)
    sorted_comp = sorted(comp)
    return sorted_comp


def max_overlap(s: str, t: str) -> int:
    """Computes the maximum overlap length between the suffix of s and the prefix of t

    Args:
        s (str): sequence
        t (str): sequence

    Returns:
        int: Max overlap length
    """
    max_l = min(len(s), len(t))
    for k in range(max_l, -1, -1):
        suffix = s[len(s)-k:]
        prefix = t[:k]
        if prefix != suffix:
            continue
        return k
    return 0


def overlapping(s: str, t: str, k: int) -> bool:
    """Returns True if there is a length k suffix of s that overlaps with a length k prefix of t, False otherwise

    Args:
        s (str): sequence
        t (str): sequence
        k (int): overlap length

    Returns:
        bool: Boolean indicating if there is overlap or not
    """
    if (len(s) < k or len(t) < k):
        return False

    suffix = s[len(s)-k:]
    prefix = t[:k]
    if suffix == prefix:
        return True
    return False


def findMotif(pattern: str, seq: str) -> dict[int, str]:
    """Finds a motif in a sequence. The motif

    Args:
        pattern (str): The motif to search, a string of the regex pattern
        seq (str): The sequence

    Returns:
        dict[int, str]: A dictionary with the 0-indexed positions of the matches as keys and the string matches
        themselves as values
    """
    motif = re.compile(f"(?=({pattern}))")
    matches: dict[int, str] = dict()
    for m in motif.finditer(seq):
        matches[m.start()] = m.group()
    return matches


def findClumps(string: str, k: int, L: int, t: int):
    """Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger) string Genome if there is an interval of Genome of length L in which Pattern appears at least t times. For example, TGCA forms a (25,3)-clump in the following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

    Args:
        string (str): _description_
        k (int): k-mer length
        L (int): maximum length of the substring containing the clumps
        t (int): minimum number of kmers needed to form clump

    Returns:
        _type_: _description_
    """
    kmers_pos_dict = defaultdict(list)
    for i, kmer in enumerate(kmer_gen(string, k)):
        kmers_pos_dict[kmer].append(i)

    kmers_pos_dict: dict[str, list[int]] = {
        kmer: positions for kmer, positions in kmers_pos_dict.items() if len(positions) >= t}
    
    clumps: set[str] = set()
    for kmer, positions in kmers_pos_dict.items():
        for i in range(len(positions)-t+1):
            L2 = positions[i+t-1] - positions[i] + k
            if L2 <= L:
                clumps.add(kmer)
    
    return clumps


def hamming_distance(seq1: str, seq2: str) -> int:
    """Computes hamming distance between two strings

    Args:
        seq1 (str): sequence1
        seq2 (str): sequence2

    Returns:
        int: _description_
    """
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def countApproximateMatches(text: str, pattern: str, d: int):
    """Given strings Text and Pattern as well as an integer d, we define Count_d(Text, Pattern) as the total number of occurrences of Pattern in Text with at most d mismatches.

    Args:
        text (str):
        pattern (str):
        d (int): maximum number of mismatches

    Returns:
        count (int): matches count
    """
    k = len(pattern)
    count: int = 0
    for i, kmer in enumerate(kmer_gen(text, k)):
        dist = hamming_distance(kmer, pattern)
        if dist <= d:
            count += 1
    return count


def approximatePatternMatching(string: str, pattern: str,  d: int):
    """We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d.

    Args:
        string (str): _description_
        pattern (str): _description_
        d (int): _description_

    Returns:
        _type_: _description_
    """
    k = len(pattern)
    idxs: list[int] = []
    for i, kmer in enumerate(kmer_gen(string, k)):
        dist = hamming_distance(kmer, pattern)
        if dist <= d:
            idxs.append(i)
    return idxs


def generate_d_neighborhood(kmer: str, d: int, alphabet: set[str]):
    """Generate all kmers that differ from kmer by at most d symbols

    Args:
        kmer (str): _description_
        d (int): _description_
        alphabet (list[str]): _description_

    Yields:
        _type_: _description_
    """
    for d1 in range(0, d+1):
        for idxs in combinations(range(len(kmer)), d1):
            alphabets = [alphabet.difference({kmer[idx]}) for idx in idxs]
            for symbs in product(*alphabets):
                kmer2_temp = [a for a in kmer]
                for idx, symb in zip(idxs, symbs):
                    kmer2_temp[idx] = symb
                kmer2 = ''.join(kmer2_temp)
                yield kmer2


def mostFrequentKmersWithMismatches(string: str, k: int, d: int, alphabet: set[str]):
    """Given strings Text and Pattern as well as an integer d, we define Countd(Text, Pattern) as the total number of occurrences of Pattern in Text with at most d mismatches. For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA, AAACA, and AAAGA. Note that two of these occurrences overlap.

    A most frequent k-mer with up to d mismatches in Text is simply a string Pattern maximizing Countd(Text, Pattern) among all k-mers. Note that Pattern does not need to actually appear as a substring of Text; for example, AAAAA is the most frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG, even though AAAAA does not appear exactly in this string.

    Args:
        string (str): _description_
        k (int): _description_
        d (int): _description_
        alphabet (set[str]): _description_

    Returns:
        _type_: _description_
    """
    kmer_freq_dict = kmers_frequency_dictionary(string, k)
    kmer_freq_with_mismatches: dict[str, int] = dict()
    for kmer, count in kmer_freq_dict.items():
        for kmer2 in generate_d_neighborhood(kmer, d, alphabet):
            if kmer2 in kmer_freq_with_mismatches:
                continue
            freq = 0
            for kmer3, count2 in kmer_freq_dict.items():
                dist = hamming_distance(kmer2, kmer3)
                if dist <= d:
                    freq += count2
            kmer_freq_with_mismatches[kmer2] = freq
    max_freq = max(kmer_freq_with_mismatches.values())
    kmers = [kmer for kmer, freq in kmer_freq_with_mismatches.items()
             if freq == max_freq]
    return kmers


def motifEnumeration(sequences: list[str], k: int, d: int, alphabet: set[str]):
    """Given a collection of strings and an integer d, a k-mer is a (k,d)-motif if it appears in every string from Dna with at most d mismatches.

    Args:
        sequences (list[str]): _description_
        k (int): _description_
        d (int): _description_
    """
    def in_seq(seq: str, kmer: str, d: int):
        for kmer2 in kmer_gen(seq, len(kmer)):
            if hamming_distance(kmer, kmer2) <= d:
                return True
        return False

    patterns: set[str] = set()
    for i, seq1 in enumerate(sequences):
        for kmer1 in get_kmers(seq1, k):
            for kmer2 in generate_d_neighborhood(kmer1, d, alphabet):
                if kmer2 in patterns:
                    continue

                in_all = all(in_seq(seq2, kmer2, d)
                             for j, seq2 in enumerate(sequences) if j != i)
                if in_all:
                    patterns.add(kmer2)
    return patterns


def getMinimumHammingDistanceToKmer(sequence: str, kmer: str) -> int:
    """Given a sequence and a kmer returns the minimum hamming distance to kmer over all k-mers of sequence

    Args:
        sequence (str): 
        kmer (str): 

    Returns:
        int: minimum hamming distance
    """
    min_d = min(hamming_distance(kmer, kmer2)
                for kmer2 in kmer_gen(sequence, len(kmer)))
    return min_d


def getMedianStrings(sequences: list[str], k: int, alphabet: set[str]):
    """Find a k-mer Pattern that minimizes d(Pattern, sequences) over all k-mers Pattern, and all sequences

    Args:
        sequences (list[str]): _description_
        k (int): _description_
        alphabet (set[str]): _description_

    Returns:
        _type_: _description_
    """
    medianStrings: set[str] = set()
    min_d = math.inf

    for kmer in generate_kmers_from_alphabet(k, sorted(alphabet)):
        dist = sum(getMinimumHammingDistanceToKmer(seq, kmer)
                   for seq in sequences)
        if dist < min_d:
            medianStrings = set([kmer])
            min_d = dist
        elif dist == min_d:
            medianStrings.add(kmer)
    return medianStrings, min_d


def generateSuffixes(text: str):
    for i in range(len(text)):
        yield text[i:]


def generateStringCycles(string: str):
    for i in range(len(string)):
        yield string[i:] + string[:i]


def suffixArray(text: str):
    suffixes = [(i, suffix) for i, suffix in enumerate(generateSuffixes(text))]
    suffixes.sort(key=lambda x: x[1])

    array = [s[0] for s in suffixes]
    return array


def partialSuffixArray(string: str, k: int) -> list[tuple[int, int]]:
    """To construct the partial suffix array SuffixArrayk(Text), we first need to construct the full suffix array and then retain only the elements of this array that are divisible by K, along with their indices i.

    Args:
        string (str): 
        k (int): 

    Returns:
        list[tuple[int, int]]: 
    """
    suffix_array = suffixArray(string)
    suffix_array = [(i, v) for i, v in enumerate(suffix_array) if v % k == 0]
    return suffix_array
