import networkx as nx
from BioInfoToolkit.Phylogeny.Phylo import PhyloTree, descendente_genotype_dist
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/mend/

A rooted binary tree can be used to model the pedigree of an individual. In this case, rather than time progressing from the root to the leaves, the tree is viewed upside down with time progressing from an individual's ancestors (at the leaves) to the individual (at the root).

An example of a pedigree for a single factor in which only the genotypes of ancestors are given is shown in Figure 1.

    Given: A rooted binary tree T in Newick format encoding an individual's pedigree for a Mendelian factor whose alleles are A (dominant) and a (recessive).

    Return: Three numbers between 0 and 1, corresponding to the respective probabilities that the individual at the root of T will exhibit the "AA", "Aa" and "aa" genotypes.
"""


def genotype_prob_from_pedigree(phylo: PhyloTree):
    tree = phylo.tree

    if not isinstance(tree, nx.DiGraph):
        raise Exception("Tree must be a digraph")

    def recurse(nodeId: int):
        nonlocal tree
        degree = tree.out_degree[nodeId]
        if degree != 0 and degree != 2:
            raise Exception("Each node must have degree 2 or 0")
        # leaf
        if degree == 0:
            genotype: str = tree.nodes[nodeId].get('name', '')
            dist: dict[str, float] = {genotype: 1.0}
            return dist
        else:
            a1, a2 = list(tree.adj[nodeId].keys())
            a1_dist = recurse(a1)
            a2_dist = recurse(a2)
            return descendente_genotype_dist(a1_dist, a2_dist)

    return recurse(0)


def verify(result: list[float], solution: list[float]) -> bool:
    correct = (len(result) == len(solution)) and all(
        abs(v1-v2) <= 0.001 for v1, v2 in zip(result, solution))
    return correct


def solve(newick: str) -> list[float]:
    phylo = PhyloTree(newick, True)

    prob_dict = genotype_prob_from_pedigree(phylo)
    sorted_keys = sorted(prob_dict.keys())

    probs = [prob_dict[key] for key in sorted_keys]
    return probs


def load_results(path: str) -> list[float]:
    lines = readTextFile(path)
    probs = [float(v) for v in lines[0].split()]
    return probs


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    newick = lines[0]

    probs = solve(newick)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(probs, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_mend_1.txt'

    lines = readTextFile(path)
    newick = lines[0]

    phylo = PhyloTree(newick, True)
    # dot = phylo.draw_dot('BT')
    # dot.view()

    prob_dict = genotype_prob_from_pedigree(phylo)
    sorted_keys = sorted(prob_dict.keys())

    probs = [prob_dict[key] for key in sorted_keys]
    print(prob_dict)

    out = ' '.join(str(prob) for prob in probs)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
