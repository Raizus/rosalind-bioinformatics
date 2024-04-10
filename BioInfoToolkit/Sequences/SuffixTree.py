from graphviz import Digraph
import networkx as nx

from BioInfoToolkit.Sequences.StringUtils import generateSuffixes

def suffix_tree_draw_dot(graph: nx.DiGraph, seq: str):
    dot = Digraph(format='svg', graph_attr={'rankdir': 'TB'})
    for nodeId in graph.nodes:
        node = graph.nodes[nodeId]
        isLeaf = graph.out_degree[nodeId] == 0
        nodeText = f"{node['start']}: {seq[node['start']:]}" if isLeaf else f"{seq[node['start']:node['end']]}"
        # nodeText = f"{node['start']}" if isLeaf else f""
        dot.node(str(nodeId), nodeText)
    for n1, nbrs in graph.adj.items():
        for n2, edge in nbrs.items():
            start, end = edge['start'], edge['end']
            label = seq[start:end]
            # label = f"s = {start}; e = {end}"
            dot.edge(str(n1), str(n2), f"{label}")

    return dot

class SuffixTree:
    # https://www.cs.jhu.edu/~langmea/resources/lecture_notes/suffix_trees.pdf
    seq: str
    root: int
    graph: nx.DiGraph

    def __init__(self, seq: str) -> None:
        self.seq = seq
        graph = nx.DiGraph()
        self.graph = graph
        n = len(seq)
        # outgoing edges, maps characters to nodes
        root = 0
        self.root = root
        graph.add_node(0, start=0, end=0, out={})
        graph.add_node(1, start=0, end=n, out={})
        graph.add_edge(0, 1, char=seq[0], start=0, end=n)
        graph.nodes[0]['out'][seq[0]] = 1
        nodeCount = 2

        for i in range(1, len(seq)):
            currId = root
            j = i
            while j < len(seq):
                if seq[j] in graph.nodes[currId]['out']:
                    childId = graph.nodes[currId]['out'][seq[j]]
                    edge = graph.edges[currId, childId]
                    start = edge['start']
                    end = edge['end']
                    label = self.seq[start:end]

                    k = j + 1
                    while k-j < len(label) and seq[k] == label[k-j]:
                        k += 1

                    if k-j == len(label):
                        currId = childId
                        j = k
                    else:
                        cExist, cNew = label[k-j], seq[k]
                        midId = nodeCount
                        newId = nodeCount + 1
                        nodeCount += 2

                        graph.remove_edge(currId, childId)

                        # update curr and child node data
                        graph.nodes[currId]['out'][seq[j]] = midId
                        currStart = graph.nodes[currId]['start']
                        currEnd = graph.nodes[currId]['end']

                        graph.add_node(midId, start=i, end=k,
                                       out={cExist: childId, cNew: newId})
                        graph.add_node(newId, start=i, end=n, out={})

                        graph.add_edge(
                            currId, midId, char=seq[j], start=j, end=k)
                        graph.add_edge(midId, childId, char=cExist,
                                       start=start+k-j, end=end)
                        graph.add_edge(midId, newId, char=cNew, start=k, end=n)

                else:
                    # make new edge and node
                    graph.add_node(nodeCount, start=i, end=n, out={})
                    graph.add_edge(currId, nodeCount,
                                   char=seq[j], start=j, end=n)
                    graph.nodes[currId]['out'][seq[j]] = nodeCount
                    nodeCount += 1

    def drawDot(self):
        dot = Digraph(format='svg', graph_attr={'rankdir': 'TB'})
        graph = self.graph
        s = self.seq
        for nodeId in graph.nodes:
            node = graph.nodes[nodeId]
            isLeaf = graph.out_degree[nodeId] == 0
            nodeText = f"{nodeId}: {s[node['start']:]}" if isLeaf else f"{nodeId}: {s[node['start']:node['end']]}"
            # nodeText = f"{node['start']}" if isLeaf else f""
            dot.node(str(nodeId), nodeText)
        for n1, nbrs in graph.adj.items():
            for n2, edge in nbrs.items():
                start, end = edge['start'], edge['end']
                label = self.seq[start:end]
                # label = f"s = {start}; e = {end}"
                dot.edge(str(n1), str(n2), f"{label}")

        return dot

    def countUniqueSubstrings(self):
        g = self.graph
        count = 0
        for n1, nbrs in g.adj.items():
            for n2, edge in nbrs.items():
                start, end = edge['start'], edge['end']
                label: str = self.seq[start:end]
                label = label.rstrip('$')
                count += len(label)
        return count

    def findLongestRepeatingSubstring(self, k: int) -> list[str]:
        """Finds the longest substrings that repeats at least k times

        Args:
            k (int): minimal length of the substring

        Returns:
            list[str]: _description_
        """
        root = self.root
        graph = self.graph
        nodes: list[int] = list(graph.nodes)

        vis = {node: False for node in nodes}

        # each node corresponds to a substring and leaf_count counts how many times each substring repeats
        leaf_count = {node: 0 for node in nodes} # count how many leaves are associated with each tree node

        def _dfs(node: int):
            vis[node] = True

            for child in graph.adj[node]:
                _dfs(child)
                leaf_count[node] += leaf_count[child]

            if len(graph.adj[node]) == 0:
                leaf_count[node] = 1

        _dfs(root)

        maxLength = 0
        result: list[str] = []
        for node in nodes:
            count = leaf_count[node]
            if count < k:
                continue

            start = graph.nodes[node]['start']
            end = graph.nodes[node]['end']
            substring = self.seq[start:end]
            l = len(substring)

            if l > maxLength:
                maxLength = l
                result = [substring]
            elif l == maxLength:
                result.append(substring)

        return result

    def findMaximalRepeats(self, minLength: int) -> list[str]:
        """A maximal repeat of a string s is a repeated substring t of s having two occurrences t1 and t2 such that t1 and t2 cannot be extended by one symbol in either direction in s and still agree.

        http://www.cs.ucf.edu/courses/cap5937/fall2004/Applications%20of%20suffix%20trees.pdf

        Args:
            k (int): minimal length of the repeats

        Returns:
            list[str]: list of maximal repeats
        """
        graph = self.graph

        maximalRepeats: list[str] = []
        node_left_diverse: dict[int, bool] = dict()
        node_left_char: dict[int, str] = dict()
        # node is left diverse <=> substring associated with node (path from root to node) is a maximal repeat
        for nodeId in nx.traversal.dfs_postorder_nodes(graph, self.root):
            node = graph.nodes[nodeId]
            if graph.out_degree[nodeId] == 0: # is leaf
                start = node['start']
                leftChar = self.seq[start-1] if start > 0 else ''
                node_left_char[nodeId] = leftChar
                continue
            
            children = list(graph.adj[nodeId])
            # is any child left diverse? if so nodeId is also left diverse
            left_diverse = any(node_left_diverse.get(child, False) for child in children)
            if left_diverse:
                node_left_diverse[nodeId] = True
                start, end = node['start'], node['end']
                substring = self.seq[start:end]
                if len(substring) >= minLength:
                    maximalRepeats.append(substring)
                continue

            chars = list(set(node_left_char[child] for child in children))
            if len(chars) == 1:
                node_left_diverse[nodeId] = False
                node_left_char[nodeId] = chars[0]
            else:
                node_left_diverse[nodeId] = True
                start, end = node['start'], node['end']
                substring = self.seq[start:end]
                if len(substring) >= minLength:
                    maximalRepeats.append(substring)

        return maximalRepeats

    def findLeafNodes(self) -> list[int]:
        graph = self.graph
        leaves: list[int] = []
        for nodeId in graph.nodes:
            od = graph.out_degree[nodeId]
            if od == 0:
                leaves.append(nodeId)
        return leaves

class SuffixTreeUkkonen:
    # https://www.cs.helsinki.fi/u/ukkonen/SuffixT1withFigs.pdf
    seq: str
    graph: nx.MultiDiGraph
    rootId: int
    nodeCount: int
    INF: int

    def addNode(self):
        nodeId = self.nodeCount
        self.graph.add_node(nodeId, out={})
        self.nodeCount += 1

        return nodeId

    def remove_edge(self, sourceId: int, destId: int, key: str):
        self.graph.remove_edge(sourceId, destId, key)
        # outDict = self.graph.nodes[sourceId]["out"]

    def get_edge_with_char(self, sNodeId: int, char: str) -> tuple[int, int, int]:
        childNodeId: int = self.graph.nodes[sNodeId]['out'][char]
        k_prime: int = self.graph.adj[sNodeId][childNodeId][char]['k']
        p_prime: int = self.graph.adj[sNodeId][childNodeId][char]['p']

        return (childNodeId, k_prime, p_prime)

    def __init__(self, t: str) -> None:
        self.INF = 10**10
        self.nodeCount = 0
        self.rootId = 0
        self.graph = nx.MultiDiGraph()
        graph = self.graph
        botId = 1
        graph.add_node(0, out={}, sLink=botId)
        graph.add_node(botId, out={})

        for i, v in enumerate(set(t)):
            graph.add_edge(botId, 0, key=v, k=-i-1, p=-i-1)
            graph.nodes[botId]['out'][v] = self.rootId

        self.nodeCount = 2
        self.seq = ' ' + t  # python indexes are 0-based, but the algorithm is 1-based
        k = 1
        sNodeId = self.rootId
        for i in range(1, len(self.seq)):
            (sNodeId, k) = self.update(sNodeId, k, i)
            (sNodeId, k) = self.canonize(sNodeId, k, i)

        graph.remove_node(1)

    def update(self, sNodeId: int, k: int, i: int):
        oldr = self.rootId
        (end_point, rNodeId) = self.test_and_split(
            sNodeId, k, i-1, self.seq[i])
        while not end_point:
            rPrimeNodeId = self.addNode()

            self.graph.add_edge(rNodeId, rPrimeNodeId,
                                key=self.seq[i], k=i, p=self.INF)
            self.graph.nodes[rNodeId]['out'][self.seq[i]] = rPrimeNodeId

            if oldr != self.rootId:
                self.graph.nodes[oldr]['sLink'] = rNodeId
            oldr = rNodeId
            sLink: int = self.graph.nodes[sNodeId]['sLink']
            (sNodeId, k) = self.canonize(sLink, k, i-1)
            (end_point, rNodeId) = self.test_and_split(
                sNodeId, k, i-1, self.seq[i])
        if oldr != self.rootId:
            self.graph.nodes[oldr]['sLink'] = sNodeId
        return (sNodeId, k)

    def test_and_split(self, sNodeId: int, k: int, p: int, char: str) -> tuple[bool, int]:
        """tests whether or not a state with canonical reference pair (s, (k, p)) is the end point, 
        that is, a state that in STrie(T i-1) would have a ti-transition.

        Args:
            nodeId (int): _description_
            k (int): _description_
            p (int): _description_
            char (str): _description_
        """
        if k <= p:
            childNodeId, k_prime, p_prime = self.get_edge_with_char(
                sNodeId, self.seq[k])
            if char == self.seq[k_prime+p-k+1]:
                return (True, sNodeId)
            else:  # make the node explicit by spliting the transition
                rNodeId = self.addNode()
                childNodeId, k_prime, p_prime = self.get_edge_with_char(
                    sNodeId, self.seq[k_prime])

                self.remove_edge(sNodeId, childNodeId, self.seq[k_prime])

                self.graph.add_edge(
                    sNodeId, rNodeId, key=self.seq[k_prime], k=k_prime, p=k_prime+p-k)
                self.graph.nodes[sNodeId]['out'][self.seq[k_prime]] = rNodeId

                self.graph.add_edge(
                    rNodeId, childNodeId, key=self.seq[k_prime + p-k+1], k=k_prime+p-k+1, p=p_prime)
                self.graph.nodes[rNodeId]['out'][self.seq[k_prime +
                                                          p-k+1]] = childNodeId

                return (False, rNodeId)
        else:
            childNodeId: int | None = self.graph.nodes[sNodeId]['out'].get(
                char, None)
            if childNodeId is None:
                return (False, sNodeId)
            return (True, sNodeId)

    def canonize(self, nodeId: int, k: int, p: int):
        if p < k:
            return (nodeId, k)
        else:
            childNodeId, k_prime, p_prime = self.get_edge_with_char(
                nodeId, self.seq[k])

            while p_prime - k_prime <= p - k:
                k += p_prime - k_prime + 1
                nodeId = childNodeId
                if k <= p:
                    childNodeId, k_prime, p_prime = self.get_edge_with_char(
                        nodeId, self.seq[k])
        return (nodeId, k)

    def draw_dot(self):
        dot = Digraph(format='svg', graph_attr={'rankdir': 'TB'})
        graph = self.graph

        for nodeId in graph.nodes:
            node = graph.nodes[nodeId]
            isLeaf = graph.out_degree[nodeId] == 0
            # nodeText = f"{node['start']}: {s[node['start']:]}" if isLeaf else f""
            dot.node(str(nodeId), f"{nodeId}")

        for n1, nbrs in graph.adj.items():
            for n2, edges in nbrs.items():
                for edgeIdx, eattr in edges.items():
                    start, end = eattr['k'], eattr['p']
                    label = ''
                    end_label = str(end) if end < self.INF else 'inf'
                    if 0 < start <= end:
                        label = self.seq[start:end+1]
                        if end == self.INF:
                            label = self.seq[start:]
                    dot.edge(str(n1), str(n2), f"{label}")

        return dot

    def countUniqueSubstrings(self):
        g = self.graph
        count = 0
        for n1, nbrs in g.adj.items():
            for n2, edges in nbrs.items():
                for edgeIdx, eattr in edges.items():
                    start, end = eattr['k'], eattr['p']
                    label: str = self.seq[start:end+1]
                    label = label.rstrip('$')
                    count += len(label)
        return count

