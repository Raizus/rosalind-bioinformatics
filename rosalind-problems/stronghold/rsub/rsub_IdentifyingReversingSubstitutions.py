from BioInfoToolkit.Phylogeny.Phylo import PhyloTree
from BioInfoToolkit.IO.IO import FASTA_lines_to_dict, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/rsub/

For a rooted tree T whose internal nodes are labeled with genetic strings, our goal is to identify reversing substitutions in T. Assuming that all the strings of T have the same length, a reversing substitution is defined formally as two parent-child string pairs (s,t) and (v,w) along with a position index i, where:
    - there is a path in T from s down to w;
    - s[i]=w[i]≠v[i]=t[i] ; and
    - if u is on the path connecting t to v, then t[i]=u[i].

In other words, the third condition demands that a reversing substitution must be contiguous: no other substitutions can appear between the initial and reversing substitution.

    Given: A rooted binary tree T with labeled nodes in Newick format, followed by a collection of at most 100 DNA strings in FASTA format whose labels correspond to the labels of T. We will assume that the DNA strings have the same length, which does not exceed 400 bp).

    Return: A list of all reversing substitutions in T (in any order), with each substitution encoded by the following three items:
        - the name of the species in which the symbol is first changed, followed by the name of the species in which it changes back to its original state
        - the position in the string at which the reversing substitution occurs; and
        - the reversing substitution in the form original_symbol->substituted_symbol->reverted_symbol.

"""

OutputT = list[tuple[str, str, int, str]]


def is_path_reversal(path: list[tuple[int, str]]) -> bool:
    if len(path) < 3:
        return False

    start_char = path[0][1]
    end_char = path[-1][1]
    if start_char == end_char:
        middle_all_equal = all(a[1] == path[1][1] for a in path[1:-1])
        if middle_all_equal and path[1][1] != start_char:
            return True
    return False


def findReversals(phylo: PhyloTree, n: int):
    tree = phylo.tree
    root = phylo.root

    reversals: list[tuple[int, list[tuple[int, str]]]] = []
    path: list[tuple[int, str]] = []

    def recurse(nodeId: int, path_index: int):
        children: list[int] = list(tree.adj[nodeId])
        char: str = tree.nodes[nodeId]['seq'][seq_index]
        path.append((nodeId, char))

        subpath = path[path_index:]
        new_path_index = path_index
        if is_path_reversal(subpath):
            reversals.append((seq_index, subpath))
            new_path_index = len(path) - 2
        elif len(subpath) == 2 and subpath[0][1] == subpath[1][1]:
            new_path_index = path_index + 1
        elif len(set(a[1] for a in subpath)) > 2:
            new_path_index = len(path) - 2

        for child in children:
            recurse(child, new_path_index)
        path.pop()

    for seq_index in range(n):
        recurse(root, 0)

    return reversals


def format_reversal(nodeid_name_dict: dict[int, str], 
                    reversal: tuple[int, list[tuple[int, str]]]) -> tuple[str, str, int, str]:
    seq_index, path = reversal
    nodeId1 = path[1][0]
    nodeId2 = path[-1][0]

    node1name = nodeid_name_dict[nodeId1]
    node2name = nodeid_name_dict[nodeId2]

    char1 = path[0][1]
    char2 = path[1][1]
    char3 = path[-1][1]

    return (node1name, node2name, seq_index+1, f"{char1}->{char2}->{char3}")


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(newick: str, fasta_dict: dict[str, str]) -> OutputT:
    phylo = PhyloTree(newick, directed=True)

    unassigned_nodes: list[int] = []
    assigned_nodes: list[int] = []

    for nodeId in phylo.tree.nodes:
        node_name = phylo.tree.nodes[nodeId]['name']
        if node_name in fasta_dict:
            seq = fasta_dict[node_name]
            phylo.tree.nodes[nodeId]['seq'] = seq
            assigned_nodes.append(nodeId)
        else:
            unassigned_nodes.append(nodeId)
            phylo.tree.nodes[nodeId]['seq'] = ''

    n = len(next(iter(fasta_dict.values())))

    # dot = phylo.draw_dot()
    # dot.view('tree1')

    reversals = findReversals(phylo, n)
    nodeid_name_dict = phylo.get_nodeId_name_dict()
    reversals = [format_reversal(nodeid_name_dict, reversal) for reversal in reversals]

    return reversals


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    reversing_subs: OutputT = []
    for line in lines:
        vals = line.split()
        n1, n2 = vals[0], vals[1]
        idx = int(vals[2])
        tran = vals[3]
        reversing_subs.append((n1, n2, idx, tran))

    return reversing_subs


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    newick = lines[0]
    fasta_data = lines[1:]
    fasta_dict = FASTA_lines_to_dict(fasta_data)

    reversals = solve(newick, fasta_dict)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(reversals, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_rsub_1.txt'

    lines = readTextFile(path)
    newick = lines[0]
    fasta_data = lines[1:]
    fasta_dict = FASTA_lines_to_dict(fasta_data)

    reversals = solve(newick, fasta_dict)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for reversal in reversals:
        out = ' '.join(str(v) for v in reversal)
        print(out)
        writeTextFile(result_path, out, 'a')

    correct = solve_and_check(path)
    print(correct)
