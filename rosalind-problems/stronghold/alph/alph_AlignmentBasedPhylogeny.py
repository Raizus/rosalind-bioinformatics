import math
from BioInfoToolkit.Phylogeny.Phylo import PhyloTree, parsimonyScore
import os
from BioInfoToolkit.IO.IO import FASTA_lines_to_dict, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/alph/

Say that we have n taxa represented by strings s1,s2,...,sn with a multiple alignment inducing corresponding augmented strings \bar{s1},\bar{s2},...,\bar{sn}.

Recall that the number of single-symbol substitutions required to transform one string into another is the Hamming distance between the strings (see “Counting Point Mutations”). Say that we have a rooted binary tree T containing \bar{s1},\bar{s2},...,\bar{sn} at its leaves and additional strings \bar{s_n+1},\bar{sn+2},...,\bar{s_2n-1} at its internal nodes, including the root (the number of internal nodes is n-1 by extension of “Counting Phylogenetic Ancestors”). Define dH(T) as the sum of dH(s¯¯¯i,s¯¯¯j) over all edges {\bar{s_i},\bar{s_j}} in T:

    dH(T)=∑{\bar{s_i},\bar{s_j}}∈E(T)dH(\bar{s_i},\bar{s_j})

Thus, our aim is to minimize dH(T).

    Given: A rooted binary tree T on n (n≤500) species, given in Newick format, followed by a multiple alignment of m (m≤n) augmented DNA strings having the same length (at most 300 bp) corresponding to the species and given in FASTA format.

    Return: The minimum possible value of dH(T), followed by a collection of DNA strings to be assigned to the internal nodes of T that will minimize dH(T) (multiple solutions will exist, but you need only output one).
"""

# the small parsimony algorithm
# https://www.coursera.org/lecture/molecular-evolution/the-small-parsimony-algorithm-Grfrh


def create_single_char_tree(phylo: PhyloTree, index: int):
    tree = phylo.tree.copy()

    for nodeId in tree.nodes:
        tree.nodes[nodeId]['seq'] = tree.nodes[nodeId]['seq'][index]
    return tree


def delta_func(a, b):
    if a == b:
        return 0
    return 1


def leafScore(target_symbol: str, symbols: list[str]):
    score: dict[str, float] = dict()

    for symb in symbols:
        score[symb] = 0 if target_symbol == symb else math.inf

    return score


def nodeScore(childrenIds: list[int], scores: dict[int, dict[str, float]], symbols: list[str]):
    parentNodeScore: dict[str, float] = dict()
    for symb in symbols:
        score = 0
        for childId in childrenIds:
            score += min([value + delta_func(symb2, symb)
                         for symb2, value in scores[childId].items()])
        parentNodeScore[symb] = score
    return parentNodeScore


def smallParsimonyAlgorithmRootedTree(phylo: PhyloTree, n: int):
    """Finds the most parsimonious (the labeling that minimizes the parsimony score) labeling of internal nodes of a rooted tree.

    Args:
        phylo (PhyloTree): _description_
        n (int): length of aligned sequences
    """
    tree = phylo.tree
    root = phylo.root
    symbols = ['A', 'C', 'G', 'T', '-']

    # maps nodeIds to scoring tables
    scores: dict[int, dict[str, float]] = dict()

    def computeScore(nodeId: int, parentNodeId: int | None, scores: dict[int, dict[str, float]], index: int):
        childrenIds: list[int] = list(tree.adj[nodeId])
        if parentNodeId is not None:
            childrenIds.remove(parentNodeId)

        if not len(childrenIds):  # leaf
            seq: str = tree.nodes[nodeId].get('seq')
            score = leafScore(seq[index], symbols)
            scores[nodeId] = score
        else:
            for childId in childrenIds:
                computeScore(childId, nodeId, scores, index)
            score = nodeScore(childrenIds, scores, symbols)
        scores[nodeId] = score

    def assignSymbol(nodeId: int, parentNodeId: int | None, scores: dict[int, dict[str, float]],
                     assignment: dict[int, str] = dict()):
        if tree.degree[nodeId] == 1:  # type: ignore
            return
        score = scores[nodeId]
        if parentNodeId is None:  # root node
            symb = ''
            min_score = math.inf
            for symb1, val in score.items():
                if val < min_score:
                    symb = symb1
                    min_score = val
            assignment[nodeId] = symb
        else:
            parentAssignment = assignment[parentNodeId]

            symb = ''
            min_score = math.inf
            for symb1, val in score.items():
                aux_score = val + delta_func(parentAssignment, symb1)
                if aux_score < min_score:
                    symb = symb1
                    min_score = aux_score
            assignment[nodeId] = symb

        childrenIds: list[int] = list(tree.adj[nodeId])
        if parentNodeId is not None:
            childrenIds.remove(parentNodeId)
        for childId in childrenIds:
            assignSymbol(childId, nodeId, scores, assignment)

    for index in range(n):
        scores: dict[int, dict[str, float]] = dict()
        computeScore(root, None, scores, index)
        # assign for each index, maps nodes to strings
        assignment: dict[int, str] = dict()
        assignSymbol(root, None, scores, assignment)
        for nodeId, symbol in assignment.items():
            tree.nodes[nodeId]['seq'] = tree.nodes[nodeId]['seq'] + symbol


InputT = tuple[str, dict[str, str]]


def verify(result: tuple[int, dict[str, str]], 
           solution: tuple[int, dict[str, str]], input: InputT) -> bool:
    dh_res = result[0]
    dh_sol = solution[0]

    # verify tree score after assignment
    # score = \sum_\{{si, sj} in E(T)\} d_H (si, sj)
    assigned_node_dict = input[1]
    unassigned_node_dict = result[1]
    phylo = PhyloTree(newick)
    for nodeId in phylo.tree.nodes:
        node_name = phylo.tree.nodes[nodeId]['name']
        if node_name in assigned_node_dict:
            seq = assigned_node_dict[node_name]
            phylo.tree.nodes[nodeId]['seq'] = seq
        else:
            seq = unassigned_node_dict[node_name]
            phylo.tree.nodes[nodeId]['seq'] = seq

    tree_score = parsimonyScore(phylo)

    correct = dh_res == dh_sol == tree_score

    return correct


def solve(newick: str, fasta_dict: dict[str, str]) -> tuple[int, dict[str, str]]:
    phylo = PhyloTree(newick)

    # dot = phylo.draw_dot()
    # dot.view('tree1')

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
    smallParsimonyAlgorithmRootedTree(phylo, n)
    score = parsimonyScore(phylo)

    unassigned_fasta_dict: dict[str,str] = dict()
    for nodeId in unassigned_nodes:
        name = phylo.tree.nodes[nodeId]['name']
        seq = phylo.tree.nodes[nodeId]['seq']
        unassigned_fasta_dict[name] = seq

    return score, unassigned_fasta_dict


def load_results(path: str) -> tuple[int, dict[str,str]]:
    lines = readTextFile(path)
    dh = int(lines[0])
    fasta_data = [line for line in lines[1:] if len(line) and not line.isspace()]
    fasta_dict = FASTA_lines_to_dict(fasta_data)
    return dh, fasta_dict


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    newick = lines[0]
    fasta_data = lines[1:]
    fasta_dict = FASTA_lines_to_dict(fasta_data)

    result = solve(newick, fasta_dict)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution, (newick, fasta_dict))
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_alph_1.txt'

    lines = readTextFile(path)
    newick = lines[0]
    fasta_data = lines[1:]
    fasta_dict = FASTA_lines_to_dict(fasta_data)

    score, unassigned_fasta_dict = solve(newick, fasta_dict)

    print(score)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(score), 'w')
    for name, seq in unassigned_fasta_dict.items():
        print(f'>{name}')
        print(seq)
        writeTextFile(result_path, f'>{name}', 'a')
        writeTextFile(result_path, seq, 'a')

    correct = solve_and_check(path)
    print(correct)

