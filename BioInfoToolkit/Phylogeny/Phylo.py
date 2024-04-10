from itertools import combinations
import math
from typing import Hashable, Literal
from graphviz import Digraph, Graph
import networkx as nx
import re

from BioInfoToolkit.Sequences.StringUtils import hamming_distance



def get_leaf_nodes(tree: nx.Graph | nx.DiGraph) -> set[int]:
    leaves = {nodeId for nodeId in tree.nodes if tree.degree[nodeId] == 1} # type: ignore
    return leaves


def descendent_dist(ancestor1: str, ancestor2: str) -> dict[str, float]:
    """Given two ancestor genotypes, returns the distribution of offspring genotypes.

    Args:
        ancestor1 (str): for example: "aa", "Aa" or "AA"
        ancestor2 (str): for example: "aa", "Aa" or "AA"

    Returns:
        dict[str, float]: offspring genotype distribution
    """
    dist: dict[str, float] = dict()
    for gene1 in ancestor1:
        for gene2 in ancestor2:
            res = ''.join(sorted(gene1 + gene2))
            dist.setdefault(res, 0.0)
            dist[res] += 1

    total = sum(dist.values())
    for key, val in dist.items():
        dist[key] = val/total
    return dist


def descendente_genotype_dist(ancestor1: dict[str, float], ancestor2: dict[str, float]) -> dict[str, float]:
    # genotypes = set(ancestor1).union(set(ancestor2))
    # genes = set(gene for genotype in genotypes for gene in genotype)
    dist: dict[str, float] = dict()

    for gt1, p1 in ancestor1.items():
        for gt2, p2 in ancestor2.items():
            res = descendent_dist(gt1, gt2)
            for g3, p3 in res.items():
                dist.setdefault(g3, 0)
                dist[g3] += p1*p2*p3

    return dist


def parse_node_str(node_str: str):
    node_str = node_str.rstrip().lstrip()
    patt = r"(\w*)(?:\s*:\s*([\d\.]+))?"
    match = re.search(patt, node_str)
    if not match:
        return '', ''
    name, weight = match.groups('')
    return name, weight


def parseNewick(newick: str, directed: bool = False):
    tree = nx.Graph()
    if directed:
        tree = nx.DiGraph()

    node_count = 0

    def createNode(name: str):
        nonlocal node_count
        node_id = node_count
        if len(name):
            tree.add_node(node_id, name=name)
        else:
            tree.add_node(node_id)

        node_count += 1
        return node_id

    node_str = ""
    node_stack: list[int] = []
    prev = ""
    for ch in reversed(newick):
        if ch == ')':
            name, weight = parse_node_str(node_str)
            nodeId = createNode(name)
            node_str = ""

            if len(node_stack):
                last_node = node_stack[-1]
                if len(weight):
                    tree.add_edge(last_node, nodeId, weight=float(weight))
                else:
                    tree.add_edge(last_node, nodeId)
            else:
                root = nodeId
            node_stack.append(nodeId)

            prev = ch

        elif ch == '(':
            if prev != ch:
                name, weight = parse_node_str(node_str)
                nodeId = createNode(name)
                node_str = ""

                last_node = node_stack[-1]
                if len(weight):
                    tree.add_edge(last_node, nodeId, weight=float(weight))
                else:
                    tree.add_edge(last_node, nodeId)
            node_stack.pop()

            prev = ch

        elif ch == ',':
            if prev != '(':
                name, weight = parse_node_str(node_str)
                nodeId = createNode(name)
                last_node = node_stack[-1]
                if len(weight):
                    tree.add_edge(last_node, nodeId, weight=float(weight))
                else:
                    tree.add_edge(last_node, nodeId)

            node_str = ""
            prev = ch

        elif ch == ';':
            node_stack = []
            prev = ch
        else:
            # ignore leading spaces
            if not (node_str == "" and ch.isspace()):
                node_str = ch + node_str
                prev = node_str

    return tree, node_count


def get_distance_to_leaves(tree: nx.Graph):
    nodes = set(nodeId for nodeId in tree.nodes)
    nnodes = len(nodes)
    leaves = get_leaf_nodes(tree)

    dist_dict: dict[int, int] = dict()

    def recurse(nodeId: int, parentNodeId: int | None):
        adj: set[int] = set(tree.adj[nodeId].keys())
        adj = adj.difference({parentNodeId})

        if len(adj) == 0:
            dist_dict[nodeId] = 0
            return 0
        else:
            dists = []
            for childNodeId in adj:
                dist = recurse(childNodeId, nodeId)
                dists.append(dist)
            min_dist = min(dists)
            dist_dict[nodeId] = min_dist
            return min_dist

    return dist_dict


def draw_phylo_dot(tree: nx.Graph | nx.DiGraph, rankdir: Literal['TB', 'BT', 'LR', 'RL'] = 'TB'):
    dot = Graph(format='svg', graph_attr={'rankdir': rankdir})
    if isinstance(tree, nx.DiGraph):
        dot = Digraph(format='svg', graph_attr={'rankdir': rankdir})

    for nodeId in tree.nodes:
        node = tree.nodes[nodeId]
        name = node.get('name', '')
        label = f"{nodeId}" + (f": {name}" if len(name) else "")
        for key, val in node.items():
            if key != 'name':
                label += f'\n{key}: {val}'
        dot.node(str(nodeId), label)

    for n1, nbrs in tree.adj.items():
        for n2, edge in nbrs.items():
            weight = edge.get('weight')
            if n2 > n1:
                if weight is not None:
                    dot.edge(str(n1), str(n2), str(weight))
                else:
                    dot.edge(str(n1), str(n2))

    return dot


def treeToNewick(tree: nx.Graph) -> str:

    def nodeToStr(nodeId: int) -> str:
        label = tree.nodes[nodeId].get('name', '')
        return label

    def recurse(nodeId: int, fromNode: int | None):
        adj = tree.adj[nodeId]
        adj = [node for node in adj if node != fromNode]

        results: list[str] = []
        for nodeId2 in adj:
            degree = tree.degree[nodeId2]  # type: ignore
            if degree == 1:
                label = nodeToStr(nodeId2)
                results.append(label)
            else:
                results.append(recurse(nodeId2, nodeId))

        return f"({','.join(results)})"

    return recurse(0, None) + ';'


def get_nodes_from_split(tree: nx.Graph | nx.DiGraph, edge: tuple[int, int]) -> tuple[set[int], set[int]]:
    """Returns the two complementary sets of nodes, resulting from cuting a tree at the given edge

    Args:
        edge (tuple[int, int]):

    Returns:
        tuple[set[int], set[int]]: S, Sc, sets of nodes in each resulting tree
    """
    tree.remove_edge(edge[0], edge[1])
    all_nodes: set[int] = set(node for node in tree.nodes)
    S: set[int] = set(nx.traversal.dfs_preorder_nodes(tree, edge[0]))
    Sc = all_nodes.difference(S)

    tree.add_edge(edge[0], edge[1])

    return S, Sc


class PhyloTree:
    tree: nx.Graph | nx.DiGraph
    node_count: int
    root = 0

    def __init__(self, newick: str, directed: bool = False) -> None:
        self.tree, self.node_count = parseNewick(newick, directed)

    def draw_dot(self, rankdir: Literal['TB', 'BT', 'LR', 'RL'] = 'TB'):
        tree = self.tree
        return draw_phylo_dot(tree, rankdir)

    def findNodeByName(self, name: str) -> int | None:
        tree = self.tree
        for nodeId in tree.nodes:
            if tree.nodes[nodeId].get('name', None) == name:
                return nodeId
        return None

    def find_non_trivial_splits(self) -> list[tuple[int, int]]:
        """Find a split S | Sc for which each set contains at least two elements. A split is given by an edge of the tree to be removed that results in a nontrivial split.

        Returns:
            list[tuple[int,int]]: list of edges (the size is the number of nontrivial splits)
        """
        tree = self.tree

        edges: list[tuple[int, int]] = []
        for n1, n2 in tree.edges:
            d1 = tree.degree[n1]  # type: ignore
            d2 = tree.degree[n2]  # type: ignore
            if d1 > 1 and d2 > 1:
                edges.append((n1, n2))

        return edges

    def get_taxon_node_dict(self) -> dict[str, int]:
        """
        Returns:
            dict[str, int]: Dictionary mapping taxon to nodeId's
        """
        tree = self.tree

        taxa = {tree.nodes[nodeId].get('name'): nodeId
                for nodeId in tree.nodes if tree.nodes[nodeId].get('name') is not None}
        return taxa

    def get_leaf_nodes(self) -> set[int]:
        tree = self.tree
        leaves = {nodeId for nodeId in tree.nodes if tree.degree[nodeId] == 1} # type: ignore
        return leaves

    def get_name_nodeId_dict(self) -> dict[str, int]:
        tree = self.tree
        name_dict = {tree.nodes[nodeId].get(
            'name', ''): nodeId for nodeId in tree.nodes if tree.nodes[nodeId].get('name') is not None}
        return name_dict

    def get_nodeId_name_dict(self) -> dict[int, str]:
        tree = self.tree
        name_dict = {nodeId: tree.nodes[nodeId]['name']
                     for nodeId in tree.nodes if tree.nodes[nodeId].get('name') is not None}
        return name_dict

    def get_leaf_dict(self) -> dict[int, str]:
        tree = self.tree
        # type: ignore
        leaves = {nodeId: tree.nodes[nodeId].get(
            'name', '') for nodeId in tree.nodes if tree.degree[nodeId] == 1}  # type: ignore
        return leaves

    def get_nodes_from_split(self, edge: tuple[int, int]) -> tuple[set[int], set[int]]:
        """Returns the two complementary sets of nodes, resulting from cuting a tree at the given edge

        Args:
            edge (tuple[int, int]):

        Returns:
            tuple[set[int], set[int]]: S, Sc, sets of nodes in each resulting tree
        """
        tree = self.tree

        tree.remove_edge(edge[0], edge[1])
        all_nodes: set[int] = set(node for node in tree.nodes)
        S: set[int] = set(nx.traversal.dfs_preorder_nodes(tree, edge[0]))
        Sc = all_nodes.difference(S)

        tree.add_edge(edge[0], edge[1])

        return S, Sc

    def get_quartets(self):
        quartets: set[tuple[tuple[int, int], tuple[int, int]]] = set()
        edges = self.find_non_trivial_splits()
        leaves = self.get_leaf_nodes()

        print("Num edges: ", len(edges))

        for i, edge in enumerate(edges):
            if i % 10 == 0:
                print("Edge: ", i)
                print("Count: ", len(quartets))
            S, Sc = self.get_nodes_from_split(edge)
            S, Sc = S.intersection(leaves), Sc.intersection(leaves)
            for pair1 in combinations(sorted(S), 2):
                for pair2 in combinations(sorted(Sc), 2):
                    if (pair1, pair2) not in quartets and (pair2, pair1) not in quartets:
                        quartets.add((pair1, pair2))

        return quartets

    def character_table_from_nontrivial_splits(self, ordered_taxa: list[str] | None = None) -> list[str]:
        character_table: list[str] = []
        edges = self.find_non_trivial_splits()

        taxa_node_dict = self.get_taxon_node_dict()
        if ordered_taxa is None:
            ordered_taxa = sorted(taxa_node_dict.keys())

        for edge in edges:
            S, Sc = self.get_nodes_from_split(edge)
            characters = [int(taxa_node_dict[taxon] not in S)
                          for taxon in ordered_taxa]

            ch_str = ''.join([str(c) for c in characters])
            character_table.append(ch_str)
        return character_table


def conflicting_splits(S1: set[int], S1c: set[int], S2: set[int], S2c: set[int]) -> bool:
    if not len(S1.intersection(S2)):
        return False
    if not len(S1.intersection(S2c)):
        return False
    if not len(S1c.intersection(S2)):
        return False
    if not len(S1c.intersection(S2c)):
        return False
    return True


def find_quartets(character_table: list[str], sort: bool = False) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    """Creates a list of quartets from a partial or full character table

    Args:
        character_table (list[str]):
        sort (bool, optional): Defaults to False.

    Returns:
        list[tuple[tuple[int, int], tuple[int, int]]]: list of quartets
    """
    quartets: set[tuple[tuple[int, int], tuple[int, int]]] = set()

    # need to check for repeats
    for character in character_table:
        idx0 = [i for i, c in enumerate(character) if c == '0']
        idx1 = [i for i, c in enumerate(character) if c == '1']

        count0 = len(idx0)
        count1 = len(idx1)

        if not (count0 >= 2 and count1 >= 2):
            continue

        for idx0_1, idx0_2 in combinations(idx0, 2):
            pair0 = (idx0_1, idx0_2)
            for idx1_1, idx1_2 in combinations(idx1, 2):
                pair1 = (idx1_1, idx1_2)
                if (pair0, pair1) not in quartets and (pair1, pair0) not in quartets:
                    quartets.add((pair0, pair1))

    quartets_list = list(quartets)
    if sort:
        quartets_list.sort(key=lambda x: x[1][1])
        quartets_list.sort(key=lambda x: x[1][0])
        quartets_list.sort(key=lambda x: x[0][1])
        quartets_list.sort(key=lambda x: x[0][0])
    return quartets_list


def newick_from_character_table(taxa: list[str], character_table: list[str]) -> str:
    """Returns an unrooted binary tree in newick format, given a taxa and consiste character table. A consistent character table is one whose characters' splits do not conflict with the edge splits of some unrooted binary tree T on the n taxa. More precisely, S1|Sc1 conflicts with S2|Sc2 if all four intersections S1 ∩ S2, S1 ∩ Sc2, Sc1 ∩ S2, and Sc1 ∩ Sc2 are nonempty.

    Args:
        taxa (list[str]): list of species / taxon
        character_table (list[str]): A consistent character table, representing the binary tree


    Returns:
        str: Unrooted binary tree in newick format
    """
    taxa = taxa.copy()
    character_table = character_table.copy()

    reduced = True
    while reduced:
        reduced = False
        for ch in character_table:
            for value in "01":
                if ch.count(value) == 2:
                    i1 = ch.find(value)
                    i2 = ch.find(value, i1 + 1)
                    taxa[i1] = f"({taxa[i1]}, {taxa[i2]})"
                    taxa = taxa[:i2] + taxa[i2+1:]
                    character_table.remove(ch)
                    reduced = True
                    for i, c in enumerate(character_table):
                        character_table[i] = character_table[i][:i2] + \
                            character_table[i][i2+1:]
                    break

    result = '(' + ','.join(taxa) + ');'
    return result


def swap_character_labels(charater: str):
    trans = charater.maketrans('10', '01')
    return charater.translate(trans)


def count_shared_splits(taxa: list[str], phylo1: PhyloTree, phylo2: PhyloTree):
    character_table1 = phylo1.character_table_from_nontrivial_splits(taxa)
    character_table2 = phylo2.character_table_from_nontrivial_splits(taxa)

    count = 0
    for ch1 in character_table1:
        for ch2 in character_table2:
            if ch1 == ch2 or ch1 == swap_character_labels(ch2):
                count += 1
                character_table2.remove(ch2)
                break
    return count


def split_distance(taxa: list[str], phylo1: PhyloTree, phylo2: PhyloTree) -> int:
    n = len(taxa)

    shared_splits = count_shared_splits(taxa, phylo1, phylo2)

    dist = 2*(n-3)-2*shared_splits

    return dist


def parsimonyScore(phylo: PhyloTree):
    """The parsimony score of the tree is the sum of hamming distances along each edge

    Args:
        phylo (PhyloTree): _description_
    """
    distance_sum = 0

    for (n1, n2) in phylo.tree.edges:
        seq1: str | None = phylo.tree.nodes[n1].get('seq', None)
        seq2: str | None = phylo.tree.nodes[n2].get('seq', None)

        if seq1 is None:
            raise Exception(
                f"Node {n1} must have a sequence associated to compute the hamming distance.")
        if seq2 is None:
            raise Exception(
                f"Node {n2} must have a sequence associated to compute the hamming distance.")

        dist = hamming_distance(seq1, seq2)
        distance_sum += dist

    return distance_sum


def count_quartets(nleaves: int):
    return math.comb(nleaves, 4)
