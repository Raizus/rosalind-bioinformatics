
from collections import Counter, OrderedDict
from itertools import accumulate, permutations
import math
from typing import TypeVar
import drawsvg as draw
import numpy as np

import networkx as nx
from graphviz import Digraph
from BioInfoToolkit.Sequences.SequenceUtils import tuple_kmers_frequency_dictionary
import re

AminoacidMonoisotopicMass = {
    "A":   71.03711,
    "C":   103.00919,
    "D":   115.02694,
    "E":   129.04259,
    "F":   147.06841,
    "G":   57.02146,
    "H":   137.05891,
    "I":   113.08406,
    "K":   128.09496,
    "L":   113.08406,
    "M":   131.04049,
    "N":   114.04293,
    "P":   97.05276,
    "Q":   128.05858,
    "R":   156.10111,
    "S":   87.03203,
    "T":   101.04768,
    "V":   99.06841,
    "W":   186.07931,
    "Y":   163.06333}

AminoacidMonoisotopicMassIntegerPlusImaginary = {
    "A":   71,
    "C":   103,
    "D":   115,
    "E":   129,
    "F":   147,
    "G":   57,
    "H":   137,
    "I":   113,
    "K":   128,
    "L":   113,
    "M":   131,
    "N":   114,
    "P":   97,
    "Q":   128,
    "R":   156,
    "S":   87,
    "T":   101,
    "V":   99,
    "W":   186,
    "Y":   163,
    "X":   4,     # imaginary
    "Z":   5}    # imaginary

ImaginaryAminoacidMonoisotopicMassInteger = {
    "X":   4,
    "Z":   5}

AminoacidMonoisotopicMassInteger = {
    "A":   71,
    "C":   103,
    "D":   115,
    "E":   129,
    "F":   147,
    "G":   57,
    "H":   137,
    "I":   113,
    "K":   128,
    "L":   113,
    "M":   131,
    "N":   114,
    "P":   97,
    "Q":   128,
    "R":   156,
    "S":   87,
    "T":   101,
    "V":   99,
    "W":   186,
    "Y":   163}

# P-E-P-F-V-A-W-A-I-Y-D-A-I-K-C-S-K-T-H-Y-N

AminoacidMonoisotopicMassOrdered = OrderedDict({
    "X":   4.0,
    "Z":   5.0,
    "G": 57.02146,
    "A": 71.03711,
    "S": 87.03203,
    "P": 97.05276,
    "V": 99.06841,
    "T": 101.04768,
    "C": 103.00919,
    "I": 113.08406,
    "L": 113.08406,
    "N": 114.04293,
    "D": 115.02694,
    "Q": 128.05858,
    "K": 128.09496,
    "E": 129.04259,
    "M": 131.04049,
    "H": 137.05891,
    "F": 147.06841,
    "R": 156.10111,
    "Y": 163.06333,
    "W": 186.07931,
})

T = TypeVar('T', int, float)

def peptideMass(peptide: str, massDict: dict[str, T] = AminoacidMonoisotopicMass):
    """Computes the mass of a protein sequence

    Args:
        protSeq (str): Protein sequence

    Returns:
        float: mass
    """
    mass = sum(massDict[aa] for aa in peptide)
    return mass


def getPossibleAAsFromMass(mass: float | int, tol: float = 0.01, 
                           aaMassDict: dict[str, int] | dict[str, float] = AminoacidMonoisotopicMass) -> list[str]:
    """Given some mass returns a list of possible aminoacids that match that mass, 
    with a margin of relative error given by tol.

    Args:
        mass (float):
        tol (float, optional): Error margin. Defaults to 0.01 (1%).

    Returns:
        list[str]: list of possible aminoacid matches, within tolerance.
    """
    aas: list[str] = []
    for aa, mass2 in aaMassDict.items():
        diff = abs(mass2-mass)
        err = diff / mass2
        if err < tol:
            aas.append(aa)
    return aas


def getAminoacidFromMass(mass: float) -> list[str]:
    """Given a target mass, returns a list with the aminoacid with the closest mass. Multiple aminoacids might be returned if they have the same mass.

    Args:
        mass (float): aminoacid mass

    Returns:
        list[str]: List of aminoacids
    """
    possibleAas: list[str] = []
    # if the masses are ordered we can find the closest
    # we could compute the differences and select the lowest
    masses = list(AminoacidMonoisotopicMassOrdered.values())
    aas = list(AminoacidMonoisotopicMassOrdered.keys())
    massesArr = np.array(masses, dtype=float)
    massesDiff = np.absolute(massesArr - mass)
    indexes = np.argsort(massesDiff)
    minDiff = massesDiff[indexes[0]]
    possibleAas.append(aas[indexes[0]])
    for i in indexes[1:]:
        massDiff = massesDiff[i]
        if massDiff == minDiff:
            possibleAas.append(aas[i])
        else:
            break
    return possibleAas


def getAminoacidsFromSpectrum(spectrum: list[float]) -> list[list[str]]:
    """Returns a aminoacid sequence that forms a protein from the prefix spectrum.

    Args:
        spectrum (list[float]): The prefix spectrum of a weighted string (collection of all its prefix weights).

    Returns:
        list[list[str]]: List of the possible aminoacids for each mass difference of the spectrum
    """
    spectrumArr = np.array(spectrum)
    diff = np.absolute(np.diff(spectrumArr))
    prot: list[list[str]] = []
    for d in diff:
        aa = getAminoacidFromMass(d)
        prot.append(aa)
    return prot


def minkowskiDiff(ms1: Counter, ms2: Counter, precision: int = 4) -> Counter:
    """Computes the Minkowski difference between two multisets.
    The Minkowski difference of multisets S1 and S2 containing real numbers is the new multiset S1 ⊖ S2 formed by taking all possible differences s1 - s2 of an element s1 from S1 and an element s2 from S2. The Minkowski sum could be defined more concisely as S1 ⊖ S2 = s1 - s2: s1 ∈ S1, s2 ∈ S2, 

    Args:
        ms1 (Counter): multiset 1
        ms2 (Counter): multiset 2
        precision (int, optional): Floating point precision to use when computing key values (s1 - s2). Defaults to 4.

    Returns:
        Counter: S1 ⊖ S2
    """

    newSet = Counter()
    for key1, val1 in ms1.items():
        for key2, val2 in ms2.items():
            newKey = round(key1 - key2, precision)
            if newKey in newSet:
                newSet[newKey] = newSet[newKey] + val1 * val2
            else:
                newSet[newKey] = val1 * val2
    orderedKeys = sorted(newSet.keys())
    orderedMultiset = Counter()
    for key in orderedKeys:
        orderedMultiset[key] = newSet[key]
    return orderedMultiset


def multisetMaxMultiplicity(ms: Counter) -> int:
    """Returns the max multiplicity of a given multiset. For a multiset S, the multiplicity of an element x is the number of times that x occurs in the set; this multiplicity is denoted S(x). Note that every set is included in the definition of multiset.

    Args:
        ms (Counter): multiset

    Returns:
        int: multiplicity
    """
    return max(ms.values())


def spectrumConvolution(spectrum: list[int]) -> Counter[int]:
    "We define the convolution of a cyclic spectrum by taking all positive differences of masses in the spectrum."
    conv: list[int] = []
    spectrum = sorted(spectrum)
    for i1, m1 in enumerate(spectrum[:-1]):
        for _, m2 in enumerate(spectrum[i1+1:]):
            diff = m2-m1
            if diff > 0:
                conv.append(diff)

    conv2 = Counter(conv)
    return conv2


def getbIons(prot: str) -> list[float]:
    """
    For example, the (unknown) protein "PRTEIN" can be cut in five possible ways: "P" and "RTEIN"; "PR" and "TEIN"; "PRT" and "EIN"; "PRTE" and "IN"; "PRTEI" and "N". We then can measure the masses of all fragments, including the entire string. The "left" end of a protein is called its N-terminus, and the ions corresponding to the protein string's prefixes (P, PR, PRT, PRTE, PRTEI) are called b-ions.

    Args:
        prot (str): Protein string

    Returns:
        list[float]: b-ions spectrum
    """
    n = len(prot)
    spectrum: list[float] = []
    prevMass = 0
    for i in range(n-1):
        aa = prot[i]
        aaMass = AminoacidMonoisotopicMass[aa]
        currMass = prevMass + aaMass
        spectrum.append(currMass)
        prevMass = currMass
    return spectrum


def getCompleteMassSpectrumPairs(prot: str, massDict: dict[str, float] | dict[str, int], includeEnds: bool = False) -> list[tuple[float, float]]:
    """Given an aminoacid sequence, returns a list of tuples with the weight of the every prefix and corresponding suffix

    Args:
        prot (str): aminoacid sequence
        includeEnds (bool, optional): Bool indicating if the pairs with empty prefix and empty suffix should be included. Defaults to False.

    Returns:
        list[tuple[float, float]]: List with the prefix and suffix mass pairs.
    """
    n = len(prot)
    spectrum: list[tuple[float, float]] = []
    protMass = peptideMass(prot)
    prevMass = 0
    if includeEnds:
        spectrum.append((0, protMass))
    # m = n if includeEnds else n - 1
    m = n
    for i in range(m):
        aa = prot[i]
        aaMass = massDict[aa]
        currMass = prevMass + aaMass
        spectrum.append((currMass, protMass - currMass))
        prevMass = currMass
    return spectrum


def getCompleteMassSpectrum(peptide: str, massDict: dict[str, float] | dict[str, int]) -> list[float] | list[int]:
    """Returns the complete mass spectrum of a protein (list of masses of all prefixes and suffixes):
        [mass(prot[0:1]), mass(prot[1:]), 
         mass(prot[0:2]), mass(prot[2:]), ..., 
         mass(prot[0:n-1]), mass(prot[n-1:])]

    Args:
        peptide (str): peptide string

    Returns:
        list[float]: mass spectrum
    """
    spectrum = []

    spectrum: list[float] = [
        mass for pair in getCompleteMassSpectrumPairs(peptide, massDict) for mass in pair]
    return spectrum


def peptideTheoreticalSpectrum(peptide: str, massDict: dict[str, int], cyclic: bool):
    """The theoretical spectrum of a peptide, is the collection of all of the masses of its subpeptides, in addition to the mass 0 and the mass of the entire peptide.

    Args:
        peptide (str): _description_
        massDict (dict[str, float] | dict[str, int]): _description_
        cyclic: Indicates if the peptide is cyclic

    Returns:
        list[float] | list[int]: _description_
    """

    peptide_tuple = tuple(massDict[aa] for aa in peptide)
    idealSpectrum: list[int] = peptideTheoreticalSpectrumTuple(peptide_tuple, cyclic)
    return idealSpectrum


def peptideTheoreticalSpectrumTuple(peptide: tuple[int,...], cyclic: bool):
    """The theoretical spectrum of a peptide, is the collection of all of the masses of its subpeptides, in addition to the mass 0 and the mass of the entire peptide.

    Args:
        peptide (tuple[int,...]): a tuple with the masses of the individual aminoacids
        massDict (dict[str, float] | dict[str, int]): _description_
        cyclic: Indicates if the peptide is cyclic

    Returns:
        list[float] | list[int]: _description_
    """
    n = len(peptide)

    spectrum = []
    for k in range(n):
        kmers_freq = tuple_kmers_frequency_dictionary(peptide, k, cyclic)
        for kmer, freq in kmers_freq.items():
            mass = sum(kmer)
            for _ in range(freq):
                spectrum.append(mass)

    spectrum.append(sum(peptide))
    return sorted(spectrum)


def peptideSpectrumScore(peptide: str, spectrum: list[int], massDict: dict[str, int], cyclic: bool) -> int:
    """Compute the score of a peptide against a spectrum.

    Args:
        peptide (str): peptide
        spectrum (list[int]): spectrum
        massDict (dict[str, int]): aminoacid mass dictionary
        cyclic: indicates if the peptide is cyclic

    Returns:
        int: score
    """
    peptide_tuple = tuple(massDict[aa] for aa in peptide)
    score = peptideSpectrumScoreTuple(peptide_tuple, spectrum, cyclic)
    return score


def peptideSpectrumScoreTuple(peptide: tuple[int, ...], spectrum: list[int], cyclic: bool) -> int:
    """Compute the score of a cyclic peptide against a spectrum.

    Args:
        peptide (str): cyclic peptide
        spectrum (list[int]): spectrum
        massDict (dict[str, int]): aminoacid mass dictionary

    Returns:
        int: score
    """
    idealSpectrum: list[int] = peptideTheoreticalSpectrumTuple(
        peptide, cyclic)

    counter1 = Counter(idealSpectrum)
    counter2 = Counter(spectrum)

    score = 0
    for m1, count1 in counter1.items():
        if m1 in counter2:
            val = min(count1, counter2[m1])
            score += val
    return score


def countPeptides(mass: int, massDict: dict[str, int] = AminoacidMonoisotopicMassInteger) -> int:
    cache: dict[int, int] = dict()

    def _recurse(mass1: int):
        if mass1 in cache:
            return cache[mass1]
        if mass1 == 0:
            return 1
        if mass1 < 0:
            return 0

        count = 0
        for aa, m1 in massDict.items():
            count += _recurse(mass1-m1)
        
        cache[mass1] = count
        return count

    total_count = _recurse(mass)
    return total_count


def getSpectrumGraph(spectrum: list[float], tol: float = 0.01) -> nx.MultiDiGraph:
    """Builds the mass spectrum graph from a spectrum (list of masses).

        For a weighted alphabet A and a collection L of positive real numbers, the spectrum graph of L is a digraph constructed in the following way. First, create a node for every real number in L. Then, connect a pair of nodes with a directed edge (u,v) if v>u and v-u is equal to the weight (w_s) of a symbol s in A. We may then label the edge with this symbol.

        If v-u in [w_s - w_s*tol, w_s + w_s*tol] matches multiple simbols then an edge will be added for each symbol, creating a directed multigraph

    Args:
        spectrum (list[float]): collection L of positive real numbers.
        tol (float, optional): _description_. Defaults to 0.01.

    Returns:
        MultiDiGraph: Spectrum Graph
    """
    n = len(spectrum)
    spectrumGraph = nx.MultiDiGraph()
    spectrumGraph.add_nodes_from(
        [(i, {"weight": w}) for i, w in enumerate(spectrum)])
    for i, j in permutations(range(n), 2):
        wi, wj = spectrum[i], spectrum[j]
        if wi < wj:
            continue
        dw = wi - wj
        aas = getPossibleAAsFromMass(dw, tol=tol)
        for aa in aas:
            aaMass = AminoacidMonoisotopicMass[aa]
            spectrumGraph.add_edge(
                j, i, label=aa, dw=dw, aaMass=aaMass, error=abs(dw-aaMass))
    return spectrumGraph


def drawSpectrumGraphDot(g: nx.MultiGraph) -> Digraph:
    """_summary_

    Args:
        g (nx.MultiGraph): Spectrum Graph

    Returns:
        Digraph: _description_
    """
    dot = Digraph(format='svg', graph_attr={'rankdir': 'LR'})
    for node in g.nodes:
        weight = g.nodes[node]["weight"]
        dot.node(str(node), f"{node}: {weight}")
    for n1, nbrs in g.adj.items():
        for n2, edges in nbrs.items():
            for edgeIdx, eattr in edges.items():
                label = eattr['label']
                dw = eattr['dw']
                aaMass = eattr['aaMass']
                error = eattr['error']
                dot.edge(str(n1), str(n2),
                            f"{label}: dw = {dw:.2f} | aaMass = {aaMass:.2f} | error = {error}")
    return dot


def longestSeqsFromSpectrumGraph(g: nx.MultiGraph) -> list[list[str]]:
    """Given the spectrum graph of a protein, returns the longest matching aminoacid sequence

    Args:
        g (nx.MultiGraph): Protein spectrum graph

    Returns:
        list[list[str]]: list with the longest aminoacid sequence
    """
    lpath = list(nx.dag_longest_path(g))  # type: ignore
    prot_seq: list[list[str]] = []
    # TODO: create path with minimum error
    for n1, n2 in zip(lpath[:-1], lpath[1:]):
        edges = g.adj[n1][n2]
        aas: list[str] = []
        for edgeIdx, eattrs in edges.items():
            label = eattrs['label']
            aas.append(label)
        prot_seq.append(aas)
    return prot_seq


class SpectrumGraph:
    graph: nx.MultiDiGraph
    spectrum: list[float]
    source: int
    sink: int

    def __init__(self, spectrum: list[float], tol: float = 0.01) -> None:
        """Builds the mass spectrum graph from a spectrum (list of masses).

        For a weighted alphabet A and a collection L of positive real numbers, the spectrum graph of L is a digraph constructed in the following way. First, create a node for every real number in L. Then, connect a pair of nodes with a directed edge (u,v) if v>u and v-u is equal to the weight (w_s) of a symbol s in A. We may then label the edge with this symbol.

        If v-u in [w_s - w_s*tol, w_s + w_s*tol] matches multiple simbols then an edge will be added for each symbol, creating a directed multigraph

        Args:
            spectrum (list[float]): collection L of positive real numbers.
            tol (float, optional): _description_. Defaults to 0.01.

        Returns:
            MultiDiGraph: Spectrum Graph
        """
        n = len(spectrum)
        graph = nx.MultiDiGraph()

        self.graph = graph
        self.spectrum = spectrum

        for i, weight in enumerate(spectrum):
            graph.add_node(i, weight=weight)

        self.source = 0
        self.sink = len(spectrum) - 1

        for i, j in permutations(range(n), 2):
            wi, wj = spectrum[i], spectrum[j]
            if wi < wj:
                continue
            dw = wi - wj
            aas = getPossibleAAsFromMass(dw, tol=tol)
            for aa in aas:
                aaMass = AminoacidMonoisotopicMass[aa]
                graph.add_edge(
                    j, i, label=aa, dw=dw, aaMass=aaMass, error=abs(dw-aaMass))

    def draw_dot(self):
        """
        Returns:
            Digraph: dot (DiGraph)
        """
        dot = drawSpectrumGraphDot(self.graph)
        return dot

    def decode_ideal_spectrum(self):
        g = self.graph
        peptide = ''
        min_error = math.inf

        sp1 = np.array(self.spectrum)

        edge_paths = nx.all_simple_edge_paths(g, self.source, self.sink)
        for edge_path in edge_paths:
            aas = [g.edges[edge]['label'] for edge in edge_path]
            _peptide = ''.join(aas)
            spectrum = sorted(getCompleteMassSpectrum(_peptide, AminoacidMonoisotopicMass))
            if len(spectrum) != len(self.spectrum):
                continue
            sp2 = np.array(spectrum)
            dist = float(np.linalg.norm(sp2-sp1))
            if dist < min_error:
                peptide = _peptide
                min_error = dist

        return peptide


def get_prefix_masses(peptide: str, aaMassDict: dict[str, int] | dict[str, float] = AminoacidMonoisotopicMass):
    prefix_masses = list(accumulate(aaMassDict[aa] for aa in peptide.upper()))
    return prefix_masses


def get_peptide_vector(peptide: str, aaMassDict: dict[str, int] | dict[str, float] = AminoacidMonoisotopicMass) -> list[int]:
    prefix_masses = get_prefix_masses(peptide, aaMassDict)
    prefix_masses = [round(mass) for mass in prefix_masses]
    max_val = prefix_masses[-1]
    peptide_vector = [0 for _ in range(max_val)]
    for mass in prefix_masses:
        peptide_vector[mass-1] = 1
    return peptide_vector


def get_peptide_from_vector(peptide_vector: list[int]) -> str:
    prefix_masses = [0] + [i+1 for i, v in enumerate(peptide_vector) if v == 1]
    diffs = [m2-m1 for m1, m2 in zip(prefix_masses, prefix_masses[1:])]

    aas = [getAminoacidFromMass(diff) for diff in diffs]
    peptide = ''.join(aa[0] for aa in aas)
    return peptide


def expandPeptides(peptides: set[str], massDict: dict[str, int]) -> set[str]:
    newSet: set[str] = set()

    for peptide in peptides:
        for aa in massDict.keys():
            peptide2 = peptide + aa
            newSet.add(peptide2)

    return newSet


def expandPeptides2(peptides: set[tuple[int, ...]], base_masses: list[int]) -> set[tuple[int, ...]]:
    newSet: set[tuple[int, ...]] = set()

    for peptide in peptides:
        for base_mass in base_masses:
            peptide2 = peptide + (base_mass,)
            newSet.add(peptide2)

    return newSet


def cyclopeptideSequencing(spectrum: list[int], massDict: dict[str, int]):
    """Given an ideal experimental spectrum, finds a cyclic peptide whose theoretical spectrum matches the experimental spectrum.

    Args:
        spectrum (list[int]): Ideal spectrum
        massDict (dict[str, int]):

    Returns:
        str | None: peptide match or None if there's no match
    """
    peptides: set[str] = set([''])
    parent_mass = max(spectrum)

    while len(peptides):
        peptides = expandPeptides(peptides, massDict)
        peptides2 = peptides.copy()

        for peptide in peptides2:
            peptide_mass = peptideMass(peptide, massDict)
            if peptide_mass == parent_mass:
                if peptideTheoreticalSpectrum(peptide, massDict, True) == spectrum:
                    return peptide
                peptides.remove(peptide)
            elif peptide_mass not in spectrum:
                peptides.remove(peptide)

    return None


def leaderboardCyclopeptideSequencing(spectrum: list[int], N: int, base_masses: list[int]):
    """_summary_

    Args:
        spectrum (list[int]): List of masses, sorted
        N (int): _description_
        base_masses (list[int]): the base masses of the spectrum, for example, corresponding to the 20 aminoacids

    Returns:
        _type_: _description_
    """
    leaderboard: dict[tuple[int,...], int] = dict()
    leaderPeptide: tuple[int,...] = ()
    leaderScore = peptideSpectrumScoreTuple(
        leaderPeptide, spectrum, False)
    leaderboard[leaderPeptide] = leaderScore
    parent_mass = max(spectrum)

    def cutLeaderboard(leaderboard: dict[tuple[int, ...], int], N: int):
        """A cut should include anyone who is tied with the Nth-place competitor. Thus, Leaderboard should be trimmed down to the “N highest-scoring peptides including ties”, which may include more than N peptides.

        Args:
            leaderboard (dict[str, int]): _description_
            N (int): _description_
        """
        scores = sorted(leaderboard.values(), reverse=True)
        if len(scores) <= N:
            return leaderboard

        Nscore = scores[N-1]
        leaderboard = {key: val for key,
                       val in leaderboard.items() if val >= Nscore}
        return leaderboard


    while len(leaderboard) > 0:
        expanded = expandPeptides2(set(leaderboard), base_masses)
        leaderboard = dict()

        for peptide in expanded:
            peptide_mass = sum(peptide)

            if peptide_mass > parent_mass:
                continue

            pep_score = leaderboard.setdefault(
                peptide, peptideSpectrumScoreTuple(peptide, spectrum, False))

            if peptide_mass == parent_mass:
                pep_score = peptideSpectrumScoreTuple(peptide, spectrum, True)
                leaderboard[peptide] = pep_score
                if pep_score > leaderScore:
                    leaderPeptide = peptide
                    leaderScore = pep_score

        leaderboard = cutLeaderboard(leaderboard, N)
    return leaderPeptide, leaderScore


def convolutionLeaderboardCyclopeptideSequencing(spectrum: list[int], M: int, N: int):
    spectrum = sorted(spectrum)
    spectral_conv = spectrumConvolution(spectrum)

    sorted_counter = sorted(
        (item for item in spectral_conv.items()), key=lambda x: x[1], reverse=True)

    sorted_spectral_conv: list[int] = []
    for mass, count in sorted_counter:
        sorted_spectral_conv.extend([mass]*count)
    sorted_spectral_conv.sort()

    # filter elements with mass >=57 ans <=200
    lb, ub = 57, 200
    sorted_counter = [item for item in sorted_counter if lb <= item[0] <= ub]

    idx = min(len(sorted_counter), M-1)
    countM = sorted_counter[idx][1]
    # top M elements with + ties with top M
    top_m_elements = sorted(
        [mass for mass, count in sorted_counter if count >= countM])

    # apply cyclopeptide leaderboard sequencing
    leaderPeptide, leaderScore = leaderboardCyclopeptideSequencing(
        spectrum, N, top_m_elements)
    return leaderPeptide, leaderScore


class SpectralGraph:
    graph: nx.MultiDiGraph
    spectral_vector: list[float] | list[int]
    source: int
    sink: int
    aaMassDict: dict[str, int] | dict[str, float]

    def __init__(self, spectral_vector: list[float] | list[int], tol: float = 0.01, 
                 aaMassDict: dict[str, int] | dict[str, float] = AminoacidMonoisotopicMassIntegerPlusImaginary) -> None:
        n = len(spectral_vector)
        graph = nx.MultiDiGraph()

        self.graph = graph
        self.spectral_vector = spectral_vector
        self.aaMassDict = aaMassDict

        for i, weight in enumerate(spectral_vector):
            graph.add_node(i, weight=weight)

        self.source = 0
        self.sink = len(spectral_vector) - 1

        graph.add_nodes_from([(i, {"weight": w}) for i, w in enumerate(spectral_vector)])
        for i, j in permutations(range(n), 2):
            dji = j - i
            aas = getPossibleAAsFromMass(dji, tol, aaMassDict)
            if len(aas):
                aa = aas[0]
                if aa not in aaMassDict:
                    continue
                aaMass = aaMassDict[aa]
                weight = graph.nodes[j]['weight']
                graph.add_edge(
                    i, j, label=aa, dw=dji, weight=-weight, aaMass=aaMass, error=abs(dji-aaMass))

    def draw_dot(self):
        """
        Returns:
            Digraph: dot (DiGraph)
        """
        dot = drawSpectrumGraphDot(self.graph)
        return dot
    
    def peptide_sequence(self):
        """
        Find amino acid string Peptide that maximizes score(peptide, spectral_vector) for all possible peptides
        """
        g = self.graph
        l = len(self.spectral_vector)
        peptide = ''
        max_score = -math.inf

        shortest_path = nx.shortest_path(g, self.source, self.sink, 'weight', 'bellman-ford')
        aas = [g.edges[(n1,n2,0)]['label'] for n1, n2 in nx.utils.pairwise(shortest_path)]
        _peptide = ''.join(aas)
        peptide = _peptide

        # edge_paths = nx.all_simple_edge_paths(g, self.source, self.sink)
        # for i, edge_path in enumerate(edge_paths):
        #     if i % 50 == 0:
        #         print(i)
        #     aas = [g.edges[edge]['label'] for edge in edge_path]
        #     _peptide = ''.join(aas)
        #     mass = int(proteinMass(_peptide, self.aaMassDict))
        #     if mass+1 != l:
        #         continue

        #     score = sum([g.nodes[edge[1]]['weight'] for edge in edge_path])

        #     if score > max_score:
        #         peptide = _peptide
        #         max_score = score

        return peptide


def peptide_score(peptide: str, spectral_vector: list[int], massDict: dict[str, float] | dict[str, int]):
    prefix_masses = get_prefix_masses(peptide, massDict)
    prefix_masses = [round(mass) for mass in prefix_masses]
    score = int(sum(spectral_vector[i]
                for i in prefix_masses if i < len(spectral_vector)))

    return score


def peptideIdentification(spectral_vector: list[int], proteome: str,
                          massDict: dict[str, float] | dict[str, int]):
    n = len(proteome)
    spvl = len(spectral_vector)
    # print("Spectral vector length: ", spvl)
    min_mass = min(val for val in massDict.values())
    max_peptide_size = math.ceil(spvl/min_mass)
    # print("Max peptide size: ", max_peptide_size)

    scoring_dict: dict[str, float] = dict()
    max_score = -math.inf
    best_fit = ''

    for i in range(n-1):
        for j in range(i+1, min(i+1+max_peptide_size, n)):
            # for j in range(i+1, n):
            peptide = proteome[i:j]
            if peptide in scoring_dict:
                continue
            mass = peptideMass(peptide, massDict)
            if mass != spvl - 1:
                continue
            score = peptide_score(peptide, spectral_vector, massDict)
            scoring_dict[peptide] = score
            if score > max_score:
                max_score = score
                best_fit = peptide
    return best_fit, max_score


def PSMSearch(spectral_vectors: list[list[int]], proteome: str, threshold: int, massDict: dict[str, float] | dict[str, int]):
    """Implements Peptide Spectrum Match Search

    Args:
        spectral_vectors (list[list[int]]): _description_
        proteome (str): _description_
        threshold (int): _description_
        massDict (dict[str, float] | dict[str, int]): _description_

    Returns:
        _type_: _description_
    """
    PSMSet: set[str] = set()
    for spectral_vector in spectral_vectors:
        peptide, score = peptideIdentification(
            spectral_vector, proteome, massDict)

        if score >= threshold:
            PSMSet.add(peptide)

    return PSMSet


def spectralDictionarySize(spectral_vector: list[int], threshold: int, max_score: int,
                           massDict: dict[str, float] | dict[str, int]):

    cache: dict[tuple[int, int], int] = dict()

    def size_func(mass: int, threshold2: int) -> int:
        if (mass, threshold2) in cache:
            return cache[(mass, threshold2)]

        if mass == 0 and threshold2 == 0:
            return 1
        # by definition threshold2<0 should be here but answer is not accepted otherwise
        elif mass <= 0 or threshold2 < 0:
            return 0

        size: int = 0
        for aa, aa_mass in massDict.items():
            si = spectral_vector[mass]
            size += size_func(round(mass-aa_mass), threshold2-si)

        cache[(mass, threshold2)] = size
        return size

    dictSize = 0
    for t in range(threshold, max_score+1):
        dictSize += size_func(len(spectral_vector)-1, t)

    return dictSize


def spectralDictionaryProbability(spectral_vector: list[int], threshold: int, max_score: int,
                                  massDict: dict[str, float] | dict[str, int]):

    cache: dict[tuple[int, int], float] = dict()
    alphabet_size = len(massDict)

    def prob_func(mass: int, threshold2: int) -> float:
        if (mass, threshold2) in cache:
            return cache[(mass, threshold2)]

        if mass == 0 and threshold2 == 0:
            return 1
        # by definition threshold2<0 should be here but answer is not accepted otherwise
        elif mass <= 0 or threshold2 < 0:
            return 0

        prob: float = 0
        for aa, aa_mass in massDict.items():
            si = spectral_vector[mass]
            prob += 1/alphabet_size * \
                prob_func(round(mass-aa_mass), threshold2-si)

        cache[(mass, threshold2)] = prob
        return prob

    dictProb = 0
    for t in range(threshold, max_score+1):
        dictProb += prob_func(len(spectral_vector)-1, t)

    return dictProb


class DiffArray:
    prefix_masses: list[int]
    prefix_mass_map: dict[int, int]

    def __init__(self, prefix_mass_map: dict[int, int], prefix_masses: list[int]) -> None:
        self.prefix_masses = prefix_masses
        self.prefix_mass_map = prefix_mass_map

    def __call__(self, mass_pref: int) -> int:
        idx = self.prefix_mass_map[mass_pref]
        return mass_pref - self.prefix_masses[idx-1]


def predecessors_previous_layer(node: tuple[int, int, int], diff: DiffArray):
    i1, j1, k1 = node
    predecessors: list[tuple[int, int, int]] = []

    if k1 == 0:
        return predecessors

    for j2 in range(0, j1):
        node2 = (i1-diff(i1), j2, k1-1)
        predecessors.append(node2)
    return predecessors


def get_valid_predecessors(node: tuple[int, int, int], diff: DiffArray, score_dict: dict[tuple[int, int, int], int]):
    i, j, t = node
    pred_k = (i-diff(i), j-diff(i), t)
    predecessors = [pred_k] + predecessors_previous_layer((i, j, t), diff)
    predecessors = [node2 for node2 in predecessors if node2 in score_dict]

    return predecessors


def findModifiedPeptide(peptide: str, spectral_vector: list[int], num_mods: int, massDict: dict[str, float] | dict[str, int]):
    """
    Given a peptide and a spectral vector, find a modified variant of this peptide that maximizes the peptide-spectrum score among all variants of the peptides with up to k modifications.

    Many proteins are subject to a large number of post-translational modifications, important for cell signaling and metabolic regulation. Most proteins are modified and the are over 600 types of aminoacid modifications

    For example adding and removing a phospate group through protein kinase and protein phosphatase (phosphorylation).

    Args:
        peptide (str):
        spectral_vector (list[int]): spectral vector starting at mass 0
        num_mods (int): maximum number of modifications
        massDict (dict[str, float] | dict[str, int]): aminoacid mass dictionary

    Returns:
        modified_peptide: str
    """
    prefix_masses = [0] + get_prefix_masses(peptide, massDict)
    prefix_masses = [int(m) for m in prefix_masses]
    peptide_mass = prefix_masses[-1]
    J = len(spectral_vector)
    # peptide_len = len(peptide)

    prefix_mass_map = {m: i for i, m in enumerate(prefix_masses)}
    diff = DiffArray(prefix_mass_map, prefix_masses)

    score_dict: dict[tuple[int, int, int], int] = dict()

    # first layer
    score_dict[(0, 0, 0)] = 0
    for i in prefix_masses[1:]:
        pred = (i-diff(i), i-diff(i), 0)
        if i < J:
            score_dict[(i, i, 0)] = spectral_vector[i] + score_dict[pred]

    for t in range(1, num_mods+1):
        for i in prefix_masses[t:]:
            for j in range(J):
                node = (i, j, t)
                predecessors = get_valid_predecessors(node, diff, score_dict)
                if len(predecessors):
                    max_val = max(score_dict[node2] for node2 in predecessors)
                    score = spectral_vector[j] + max_val
                    score_dict[node] = score

    end_nodes = list((peptide_mass, J-1, t) for t in range(num_mods+1)
                     if (peptide_mass, J-1, t) in score_dict)
    max_val = max(score_dict[node] for node in end_nodes)

    max_nodes = [node for node in end_nodes if score_dict[node] == max_val]

    # backtrack
    node = max_nodes[0]
    path = [node]
    while node != (0, 0, 0):
        i, j, t = node
        predecessors = get_valid_predecessors(node, diff, score_dict)
        for pred in predecessors:
            if score_dict[node] == spectral_vector[j] + score_dict[pred]:
                path.insert(0, pred)
                node = pred
                break

    # construct modified peptide from path
    peptide_mod = ''
    for i, (n1, n2) in enumerate(zip(path, path[1:])):
        peptide_mod += peptide[i]
        i1, j1, k1 = n1
        i2, j2, k2 = n2
        delta = (j2-j1) - (i2-i1)
        if k1 != k2:
            if delta > 0:
                peptide_mod += f"(+{delta})"
            elif delta < 0:
                peptide_mod += f"({delta})"

    return peptide_mod


class PSMGraph:
    graph: nx.DiGraph
    peptide: str
    spectral_vector: list[int]
    prefix_masses: list[int]
    prefix_mass_map: dict[int, int]

    def __init__(self, peptide: str, spectral_vector: list[int], massDict: dict[str, float] | dict[str, int]) -> None:
        prefix_masses = [0] + get_prefix_masses(peptide, massDict)
        prefix_masses = [int(m) for m in prefix_masses]
        peptide_mass = prefix_masses[-1]
        spcvl = len(spectral_vector)
        peptide_len = len(peptide)
        prefix_mass_map = {m: i for i, m in enumerate(prefix_masses)}

        graph = nx.DiGraph()

        self.graph = graph
        self.peptide = peptide
        self.spectral_vector = spectral_vector
        self.prefix_masses = prefix_masses
        self.prefix_mass_map = prefix_mass_map

        for prefix_mass in prefix_masses:
            for i, si in enumerate(spectral_vector):
                graph.add_node((prefix_mass, i), weight=si)

        prev_nodes: list[tuple[int, int]] = [(0, 0)]
        for i, prefix_mass in enumerate(prefix_masses[:-2]):
            next_nodes: list[tuple[int, int]] = []
            i2 = self.prefix_masses[i+1]
            for node1 in prev_nodes:
                i1, j1 = node1
                for j2 in range(j1+1, spcvl-peptide_len+i+1):
                    node2 = (i2, j2)
                    graph.add_edge(node1, node2)
                    next_nodes.append(node2)
            prev_nodes = next_nodes

        i2 = self.prefix_masses[-1]
        for node1 in prev_nodes:
            j2 = spcvl - 1
            node2 = (i2, j2)
            graph.add_edge(node1, node2)

    def diff(self, mass_pref: int) -> int:
        idx = self.prefix_mass_map[mass_pref]
        return mass_pref - self.prefix_masses[idx-1]

    def draw_svg(self):
        cell_size = 50
        font_size = 15
        n_rows = len(self.peptide) + 1
        n_cols = len(self.spectral_vector)
        height = cell_size*(n_rows+3)
        width = cell_size*(n_cols+3)
        d = draw.Drawing(width, height)
        d.append(draw.Rectangle(0, 0, width, height, fill="white"))

        # draw nodes
        node_r = 0.15*cell_size
        for row in range(n_rows):
            for col in range(n_cols):
                c = draw.Circle((col + 2.5)*cell_size, (row + 2.5)
                                * cell_size, node_r, fill='lightgray')
                d.append(c)

        # draw mass prefix labels
        for i, prefix_mass in enumerate(self.prefix_masses):
            text = draw.Text(str(prefix_mass), font_size, (1.5)*cell_size,
                             (i + 2.5)*cell_size, text_anchor='middle', center=True)
            d.append(text)

        # draw aminoacid labels
        for i, aa in enumerate(self.peptide):
            text = draw.Text(aa, font_size, (0.5)*cell_size,
                             (i + 1 + 2.5)*cell_size, text_anchor='middle', center=True)
            d.append(text)

        # draw column mass labels
        for j in range(len(self.spectral_vector)):
            text = draw.Text(str(j), font_size, (j + 2.5)*cell_size,
                             (1.5)*cell_size, text_anchor='middle', center=True)
            d.append(text)

        # draw edges
        for edge in self.graph.edges:
            n1, n2 = edge
            i1 = self.prefix_mass_map[n1[0]]
            i2 = self.prefix_mass_map[n2[0]]
            p1 = np.array([(n1[1] + 2.5)*cell_size, (i1 + 2.5)*cell_size])
            p2 = np.array([(n2[1] + 2.5)*cell_size, (i2 + 2.5)*cell_size])
            p12 = p2 - p1
            p12 = p12 / np.linalg.norm(p12)
            p1 = p1 + node_r * p12
            p2 = p2 - node_r * p12

            line = draw.Line(p1[0], p1[1], p2[0], p2[1],
                             stroke_width=2, stroke='lightblue')
            d.append(line)

        return d


def parseModifiedPeptideToArray(mod_peptide: str) -> list[int]:
    masses: list[int] = []
    mass_dict = AminoacidMonoisotopicMassInteger
    for match in re.finditer(r'([A-Z])(?:\(([+-]\d+)\))?', mod_peptide):
        aa = match[1]
        mod = int(match[2]) if match.group(2) is not None else 0
        mass = mass_dict[aa] + mod
        masses.append(mass)
    return masses
