from BioInfoToolkit.Phylogeny.Phylo import PhyloTree
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/ctbl/

Given a collection of n taxa, any subset S of these taxa can be seen as encoding a character that divides the taxa into the sets S and Sc; we can represent the character by S|Sc, which is called a split. Alternately, the character can be represented by a character array A of length n for which A[j]=1 if the jth taxon belongs to S and A[j]=0 if the jth taxon belongs to Sc (recall the "ON"/"OFF" analogy from “Counting Subsets”).

At the same time, observe that the removal of an edge from an unrooted binary tree produces two separate trees, each one containing a subset of the original taxa. So each edge may also be encoded by a split S|Sc.

A trivial character isolates a single taxon into a group of its own. The corresponding split S|Sc
must be such that S or Sc contains only one element; the edge encoded by this split must be incident to a leaf of the unrooted binary tree, and the array for the character contains exactly one 0 or exactly one 1. Trivial characters are of no phylogenetic interest because they fail to provide us with information regarding the relationships of taxa to each other. All other characters are called nontrivial characters (and the associated splits are called nontrivial splits).

A character table is a matrix C in which each row represents the array notation for a nontrivial character. That is, entry Ci,j denotes the "ON"/"OFF" position of the ith character with respect to the jth taxon.

    Given: An unrooted binary tree T in Newick format for at most 200 species taxa.

    Return: A character table having the same splits as the edge splits of T. The columns of the character table should encode the taxa ordered lexicographically; the rows of the character table may be given in any order. Also, for any given character, the particular subset of taxa to which 1s are assigned is arbitrary.
"""


def verify(result: list[str], solution: list[str]) -> bool:
    correct = set(result) == set(solution)
    return correct


def solve(newick: str) -> list[str]:
    phylo = PhyloTree(newick)
    character_table = phylo.character_table_from_nontrivial_splits()
    return character_table


def load_results(path: str) -> list[str]:
    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    return lines


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    newick = lines[0]
    character_table = solve(newick)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(character_table, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ctbl_1.txt'

    lines = readTextFile(path)
    newick = lines[0]
    phylo = PhyloTree(newick)

    # dot = phylo.draw_dot()
    # dot.view()

    character_table = phylo.character_table_from_nontrivial_splits()

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for c in character_table:
        print(c)
        writeTextFile(result_path, c, 'a')

    correct = solve_and_check(path)
    print(correct)
