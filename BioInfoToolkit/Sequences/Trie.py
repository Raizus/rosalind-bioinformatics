import networkx as nx

# class Trie:
#     graph: nx.DiGraph
#     rootId: int
#     nodeCount: int

#     def __init__(self, patterns: list[str]) -> None:
#         graph = nx.DiGraph()

#         self.rootId = 0
#         graph.add_node(0)
#         self.nodeCount = 1

#         for string in patterns:
#             currentNode = self.rootId
#             for i, char in enumerate(string):
#                 pass


#     def draw_dot(self):
#         pass


from typing import Generic, Hashable, TypeVar
from graphviz import Digraph

from pyparsing import Iterable

T = TypeVar("T")


class TrieNode:
    children: dict[Hashable, "TrieNode"]
    id: int
    isEnd: bool

    def __init__(self, id: int):
        self.children = dict()
        self.isEnd = False
        self.id = id


class Trie:
    root: TrieNode
    lastId: int

    def __init__(self) -> None:
        self.root = TrieNode(0)
        self.lastId = 0

    def insert(self, key: Iterable):
        currentNode = self.root
        for item in key:
            if item not in currentNode.children:
                self.lastId += 1
                currentNode.children[item] = TrieNode(self.lastId)
            currentNode = currentNode.children[item]
        currentNode.isEnd = True

    def draw_dot(self):
        dot = Digraph(format='svg', graph_attr={'rankdir': 'TB'})

        def build(node: TrieNode):
            dot.node(str(node.id), str(node.id))
            for key, childNode in node.children.items():
                dot.edge(str(node.id), str(childNode.id), str(key))
            for childNode in node.children.values():
                build(childNode)

        build(self.root)

        return dot

    def prefix_trie_matching(self, text: str) -> str | None:
        node = self.root
        match = ''
        i = 0
        symbol = text[i]
        while True:
            if len(node.children) == 0:
                return match

            next_node = node.children.get(symbol)
            if next_node:
                match = match + symbol
                node = next_node    
                i+=1

                if len(node.children) == 0:
                    return match
                symbol = text[i]
            else: # no match
                return None
        
    def trie_matching(self, text: str):
        idxs: list[int] = []
        for i in range(len(text)):
            match = self.prefix_trie_matching(text[i:])
            if match is not None:
                idxs.append(i)
        return idxs


class TrieNx:
    graph: nx.DiGraph
    nodeCount: int
    root: int

    def __init__(self) -> None:
        graph = nx.DiGraph()
        self.graph = graph
        self.nodeCount = 0
        self.root = 0

        # root
        graph.add_node(self.nodeCount, children={})
        self.nodeCount += 1

    def insert(self, pattern: str):
        g = self.graph
        currentNode = self.root
        for symbol in pattern:
            children = g.nodes[currentNode]['children']
            # children_dict = {edge[char]: node for node, edge in succ.items()}
            if symbol not in children:
                g.add_node(self.nodeCount, children={})
                g.add_edge(currentNode, self.nodeCount, char=symbol)
                children[symbol] = self.nodeCount
                self.nodeCount += 1
            currentNode = children[symbol]
            
    def draw_dot(self):
        dot = Digraph(format='svg', graph_attr={'rankdir': 'TB'})
        g = self.graph
        for node in g.nodes:
            dot.node(str(node), str(node))
        for edge in g.edges:
            char = g.edges[edge]['char']
            dot.edge(str(edge[0]), str(edge[1]), char)            

        return dot
    
    def prefix_trie_matching(self, text: str) -> str | None:
        current_node = self.root
        graph = self.graph
        match = ''
        for i, symbol in enumerate(text):
            if len(graph.succ[current_node]) == 0:
                return match

            next_node = graph.nodes[current_node]['children'].get(symbol)
            if next_node:
                match = match + symbol
                current_node = next_node

                if len(graph.succ[current_node]) == 0:
                    return match
            else:  # no match
                return None
        return None

    def trie_matching(self, text: str):
        idxs: list[int] = []
        for i in range(len(text)):
            match = self.prefix_trie_matching(text[i:])
            if match is not None:
                idxs.append(i)
        return idxs
