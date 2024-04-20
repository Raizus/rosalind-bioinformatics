from BioInfoToolkit.Phylogeny.Phylo import PhyloTree, newick_from_character_table
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os
import networkx as nx


"""
https://rosalind.info/problems/chbp/

Because a tree having n nodes has n-1 edges (see “Completing a Tree”), removing a single edge from a tree will produce two smaller, disjoint trees. Recall from “Creating a Character Table” that for this reason, each edge of an unrooted binary tree corresponds to a split S|Sc, where S is a subset of the taxa.

A consistent character table is one whose characters' splits do not conflict with the edge splits of some unrooted binary tree T
on the n taxa. More precisely, S1|Sc1 conflicts with S2|Sc2 if all four intersections S1∩S2, S1∩Sc2, Sc1∩S2, and Sc1∩Sc2 are nonempty. As a simple example, consider the conflicting splits {a,b}|{c,d} and {a,c}|{b,d}.

More generally, given a consistent character table C, an unrooted binary tree T "models" C if the edge splits of T agree with the splits induced from the characters of C.

    Given: A list of n species (n≤80) and an n-column character table C in which the jth column denotes the jth species.

    Return: An unrooted binary tree in Newick format that models C.
"""


def verify(result: str, solution: str) -> bool:
    phylo1 = PhyloTree(result)
    phylo2 = PhyloTree(solution)

    correct = nx.is_isomorphic(phylo1.tree, phylo2.tree)
    return correct


def solve(taxa: list[str], character_table: list[str]) -> str:
    newick = newick_from_character_table(taxa, character_table)
    return newick


def load_results(path: str) -> str:
    lines = readTextFile(path)
    newick = lines[0]
    return newick


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    taxa = lines[0].split(' ')
    character_table = lines[1:]

    newick = solve(taxa, character_table)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(newick, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_chbp_1.txt'

    lines = readTextFile(path)
    taxa = lines[0].split(' ')
    character_table = lines[1:]

    newick = newick_from_character_table(taxa, character_table)
    print(newick)

    # phylo2 = PhyloTree(newick)
    # dot = phylo2.draw_dot()
    # dot.view('tree1')

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, newick, 'w')

    # correct = solve_and_check(path)
    # print(correct)
