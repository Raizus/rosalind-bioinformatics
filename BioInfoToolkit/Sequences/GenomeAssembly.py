from abc import ABC, abstractmethod
from collections import defaultdict
from functools import reduce
from itertools import accumulate, combinations
from typing import Iterable

from graphviz import Digraph
from BioInfoToolkit.Sequences.StringUtils import kmer_gen, max_overlap, overlapping
import networkx as nx

def overlap_graph_assembly(sequences: dict[str, str]) -> defaultdict[str, dict[str, int]]:
    """Builds and returns an overlap graph (directed graph). For a collection of strings, the overlap graph for the strings is a directed graph O_k in which each string is represented by a node, and string s is connected to 
    string t with a directed edge when there is a length k suffix of s that matches a
    length k prefix of t, as long as s≠t; we demand s≠t to prevent directed loops in the 
    overlap graph (although directed cycles may be present), where k > len(s)/2 and k > len(t)/2

    Args:
        sequences (dict[str, str]): dictionary of sequences, where keys are sequence ID's and values are sequences

    Returns:
        defaultdict[str, set[str]]: The overlap graph
    """
    graph: defaultdict[str, dict[str, int]] = defaultdict(dict)
    for (key1, key2) in combinations(sequences.keys(), 2):
        seq1 = sequences[key1]
        seq2 = sequences[key2]
        if seq1 == seq2:
            continue
        a = max_overlap(seq1, seq2)
        if a > len(seq1)/2 and a > len(seq2)/2:
            graph[key1][key2] = a
        b = max_overlap(seq2, seq1)
        if b > len(seq1)/2 and b > len(seq2)/2:
            graph[key2][key1] = b
    return graph


def superstring_from_overlap_graph(sequences: dict[str, str], ol_graph: defaultdict[str, dict[str, int]]):
    # find starting nodes (outgoing node)
    start_nodes = []
    for node in ol_graph.keys():
        not_in_edges = True
        for edges in ol_graph.values():
            if node in edges:
                not_in_edges = False
                break
        if not_in_edges:
            start_nodes.append(node)

    start_node = start_nodes[0]
    path = [start_node]
    node = start_node
    while node in ol_graph:
        edges = ol_graph[node]
        edge = next(iter(edges))
        path.append(edge)
        node = edge

    superstring = sequences[path[0]]
    for n1, n2 in zip(path, path[1:]):
        k = ol_graph[n1][n2]
        superstring += sequences[n2][k:]
    return superstring


class OverlapGraph:
    """Given an arbitrary collection of k-mers Patterns, we form a graph having a node for each k-mer in Patterns and connect k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern'). The resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns).
    """
    graph: nx.DiGraph
    k: int

    def __init__(self, sequences: dict[str, str] | list[str], k: int) -> None:
        """Builds and returns an overlap graph (directed graph). For a collection of strings 
        and a positive integer k, the overlap graph for the strings is a directed graph 
        O_k in which each string is represented by a node, and string s is connected to 
        string t with a directed edge when there is a length k suffix of s that matches a
        length k prefix of t, as long as s≠t; we demand s≠t to prevent directed loops in the 
        overlap graph (although directed cycles may be present).

        Args:
            sequences (dict[str, str]): dictionary of sequences, where keys are sequence ID's and values are sequences
            k (int): overlap length

        Returns:
            defaultdict[str, set[str]]: The overlap graph
        """
        graph = nx.DiGraph()
        self.graph = graph
        self.k = k

        if type(sequences) == list:
            n = len(sequences)

            for id, sequence in enumerate(sequences):
                graph.add_node(str(id), seq=sequence)

            for id1, id2 in combinations(range(n), 2):
                seq1 = sequences[id1]
                seq2 = sequences[id2]

                if seq1 == seq2:
                    continue

                overlaps = overlapping(seq1, seq2, k)
                if overlaps:
                    graph.add_edge(str(id1), str(id2))

                overlaps = overlapping(seq2, seq1, k)
                if overlaps:
                    graph.add_edge(str(id2), str(id1))

        elif type(sequences) == dict:

            for id, sequence in sequences.items():
                graph.add_node(str(id), seq=sequence)

            for (key1, key2) in combinations(sequences.keys(), 2):
                seq1 = sequences[key1]
                seq2 = sequences[key2]
                if seq1 == seq2:
                    continue

                overlaps = overlapping(seq1, seq2, k)
                if overlaps:
                    graph.add_edge(key1, key2)

                overlaps = overlapping(seq2, seq1, k)
                if overlaps:
                    graph.add_edge(key2, key1)

    def draw_dot(self):
        dot = Digraph(format='svg', graph_attr={'rankdir': 'LR'})
        g = self.graph
        for node in g.nodes:
            sequence = g.nodes[node]["seq"]
            dot.node(str(node), f"{node}: {sequence}")

        for n1, nbrs in g.adj.items():
            for n2, edge in nbrs.items():
                dot.edge(str(n1), str(n2))
        return dot


def deBruijnMultiGraphFromString(string: str, k: int):
    graph = nx.MultiDiGraph()
    for i, kmer in enumerate(kmer_gen(string, k)):
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if prefix not in graph.nodes:
            graph.add_node(prefix)
        if suffix not in graph.nodes:
            graph.add_node(suffix)

        graph.add_edge(prefix, suffix, seq=kmer)    
    return graph


def deBruijnMultiGraphFromKmers(kmers: Iterable[str]):
    graph = nx.MultiDiGraph()
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if prefix not in graph.nodes:
            graph.add_node(prefix)
        if suffix not in graph.nodes:
            graph.add_node(suffix)

        graph.add_edge(prefix, suffix, seq=kmer)    
    return graph


def deBruijnMultiGraphFromPairedKmers(kmerPairs: Iterable[tuple[str, str]]):
    graph = nx.MultiDiGraph()
    for pair in kmerPairs:
        prefix = (pair[0][:-1], pair[1][:-1])
        suffix = (pair[0][1:], pair[1][1:])

        if prefix not in graph.nodes:
            graph.add_node(prefix)
        if suffix not in graph.nodes:
            graph.add_node(suffix)

        graph.add_edge(prefix, suffix, seq=pair)
    return graph


class DeBruijnMultiGrapAbstract(ABC):
    graph: nx.MultiDiGraph
    k: int

    def draw_dot(self):
        dot = Digraph(format='svg', graph_attr={'rankdir': 'LR'})
        g = self.graph
        for node in g.nodes:
            dot.node(str(node), f"{node}")

        for n1, nbrs in g.adj.items():
            for n2, edges in nbrs.items():
                for edgeIdx, edgeData in edges.items():
                    seq = g.edges[(n1, n2, edgeIdx)]['seq']
                    dot.edge(str(n1), str(n2), f"{seq}")
        return dot

    def hasEulerianPath(self) -> bool:
        return nx.has_eulerian_path(self.graph)
    
    def hasEulerianCycle(self) -> bool:
        return nx.is_eulerian(self.graph)

    def is_balanced(self):
        g = self.graph
        return all((g.in_degree[node] == g.out_degree[node] for node in g.nodes))

    @abstractmethod
    def reconstructStringFromEulerianPath(self) -> str:
        pass

    @abstractmethod
    def reconstructStringFromEulerianCycle(self) -> str:
        pass


class DeBruijnMultiGraph(DeBruijnMultiGraphAbstract):
    def __init__(self, graph: nx.MultiDiGraph, k: int) -> None:
        self.graph = graph
        self.k = k

    def reconstructStringFromEulerianPath(self):
        edgesPathGen = nx.eulerian_path(self.graph, keys=True)

        string = ''
        firstEdge = next(edgesPathGen)
        label = self.graph.edges[firstEdge]['seq']
        string += label

        for edge in edgesPathGen:
            label = self.graph.edges[edge]['seq']
            string += label[-1]
        return string

    def reconstructStringFromEulerianCycle(self):
        edgesPathGen = nx.eulerian_circuit(self.graph, keys=True)

        string = ''
        firstEdge = next(edgesPathGen)
        label = self.graph.edges[firstEdge]['seq']
        string += label

        for edge in edgesPathGen:
            label = self.graph.edges[edge]['seq']
            string += label[-1]
        return string[:-self.k+1]

    def generateMaximalNonBranchingPaths(self):
        def in1out1node(node: str) -> bool:
            out_degree = g.out_degree[node]
            in_degree = g.in_degree[node]

            return in_degree == 1 and out_degree == 1

        g = self.graph
        unvisitedNodes = set(g.nodes)
        for node1 in g.nodes:
            out_degree = g.out_degree[node1]

            if in1out1node(node1):
                continue

            if out_degree == 0:
                continue

            adjacencyDict = g.adj[node1]
            unvisitedNodes.discard(node1)
            for node2, edgesDict in adjacencyDict.items():
                unvisitedNodes.discard(node2)
                for edgeIdx, edgesData in edgesDict.items():
                    newPath: list[tuple[str,str,int]] = [(node1, node2, edgeIdx)]
                    while in1out1node(node2):
                        node3 = next(g.successors(node2))
                        newPath.append((node2, node3, 0))
                        node2 = node3
                        unvisitedNodes.discard(node3)
                    yield newPath

        # add isolated cycles
        if len(unvisitedNodes):
            print("There are isolated cycles")
        
    def edgesPathToString(self, edges: list[tuple[str,str,int]]):
        string = ''

        for i, edge in enumerate(edges):
            label = self.graph.edges[edge]['seq']
            if i == 0:
                string += label
            else:
                string += label[-1]
        return string

    def generateContigs(self):
        for path in self.generateMaximalNonBranchingPaths():
            contig = self.edgesPathToString(path)
            yield contig

    def findAllEulerianCycles(self, start_node_id: str):
        num_edges = len(self.graph.edges)
        g = self.graph

        seqs: list[str] = []

        def recurse(node_id: str, path: list[str], count: int):
            adj = dict(g.adj[node_id])
            path.append(node_id)
            count += 1

            if len(adj) == 0:
                if count == num_edges+1 and node_id == start_node_id: # found the eulerian path
                    seq = reduce(lambda prev,next: prev+next[-1], path[1:], path[0])
                    l=len(node_id)
                    seqs.append(seq[:-l])
            else:
                for child_id, edge_dict in adj.items():
                    edge_dict = dict(edge_dict)
                    g.remove_edge(node_id, child_id)
                    recurse(child_id, path, count)
                    g.add_edge(node_id, child_id, seq=node_id+child_id[-1])

            path.pop()

        recurse(start_node_id, [], 0)
        a = 0


class PairedDeBruijnMultiGraph(DeBruijnMultiGraphAbstract):
    d: int

    def __init__(self, graph: nx.MultiDiGraph, k: int, d: int) -> None:
        self.graph = graph
        self.k = k
        self.d = d
    
    def reconstructStringFromEulerianPath(self):
        edgesPathGen = nx.eulerian_path(self.graph, keys=True)

        string1 = ''
        string2 = ''
        firstEdge = next(edgesPathGen)
        label = self.graph.edges[firstEdge]['seq']
        string1 += label[0]
        string2 += label[1]

        for edge in edgesPathGen:
            label = self.graph.edges[edge]['seq']
            string1 += label[0][-1]
            string2 += label[1][-1]
        
        d, k = self.d, self.k
        string = string1 + string2[len(string2)-(d):]
        return string

    def reconstructStringFromEulerianCycle(self):
        raise NotImplementedError()
