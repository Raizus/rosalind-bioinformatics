
from itertools import combinations
import math
from BioInfoToolkit.Phylogeny.Phylo import PhyloTree, count_quartets, count_shared_splits, swap_character_labels
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/qrtd/

In “Counting Quartets”, we found an expression for q(T), the number of quartets that can be inferred from an unrooted binary tree containing n taxa.

If T1 and T2 are both unrooted binary trees on the same n taxa, then we now let q(T1,T2) denote the number of inferred quartets that are common to both trees. The quartet distance between T1 and T2, dq(T1,T2) is the number of quartets that are only inferred from one of the trees. More precisely, dq(T1,T2)=q(T1)+q(T2)-2q(T1,T2).

    Given: A list containing n taxa (n≤2000) and two unrooted binary trees T1 and T2 on the given taxa. Both T1 and T2 are given in Newick format.

    Return: The quartet distance dq(T1,T2).
"""


def count_shared_quartets2(taxa: list[str], phylo1: PhyloTree, phylo2: PhyloTree):
    char_table_1 = set(phylo1.character_table_from_nontrivial_splits(taxa))
    char_table_2 = set(phylo2.character_table_from_nontrivial_splits(taxa))

    nleaves = len(taxa)
    quartets: set[tuple[tuple[int, int], tuple[int, int]]] = set()

    ref_set = set(range(nleaves))

    def find_shared(ch1: str, ch2: str):
        F1 = {i for i, c in enumerate(ch1) if c == '1'}
        F1c = ref_set.difference(F1)

        G1 = {i for i, c in enumerate(ch2) if c == '1'}
        G1c = ref_set.difference(G1)

        A1 = F1.intersection(G1)
        A1c = F1c.intersection(G1c)

        B1 = F1.intersection(G1c)
        B1c = F1c.intersection(G1)

        shared: set[tuple[tuple[int, int], tuple[int, int]]] = set()

        for comb1 in combinations(A1, 2):
            pair1 = comb1 if comb1[0] < comb1[1] else (comb1[1], comb1[0])
            for comb2 in combinations(A1c, 2):
                pair2 = comb2 if comb2[0] < comb2[1] else (comb2[1], comb2[0])
                if (pair1, pair2) not in quartets and (pair1, pair2) not in quartets:
                    quartets.add((pair1, pair2))

        for comb1 in combinations(B1, 2):
            pair1 = comb1 if comb1[0] < comb1[1] else (comb1[1], comb1[0])
            for comb2 in combinations(B1c, 2):
                pair2 = comb2 if comb2[0] < comb2[1] else (comb2[1], comb2[0])
                if (pair1, pair2) not in quartets and (pair1, pair2) not in quartets:
                    quartets.add((pair1, pair2))

    for i, ch1 in enumerate(char_table_1):
        print(i)
        for j, ch2 in enumerate(char_table_2):
            find_shared(ch1, ch2)
    return len(quartets)


def count_shared_quartets(taxa: list[str], phylo1: PhyloTree, phylo2: PhyloTree):
    char_table_1 = set(phylo1.character_table_from_nontrivial_splits(taxa))
    char_table_2 = set(phylo2.character_table_from_nontrivial_splits(taxa))

    nleaves = len(taxa)
    ref_set = set(range(nleaves))

    def shared(ch1: str, ch2: str):
        F1 = {i for i, c in enumerate(ch1) if c == '1'}
        F1c = ref_set.difference(F1)

        G1 = {i for i, c in enumerate(ch2) if c == '1'}
        G1c = ref_set.difference(G1)

        A = math.comb(len(F1.intersection(G1)), 2)
        B = math.comb(len(F1.intersection(G1c)), 2)
        C = math.comb(len(F1c.intersection(G1)), 2)
        D = math.comb(len(F1c.intersection(G1c)), 2)

        return A*D+B*C

    count = 0
    for ch1 in char_table_1:
        for ch2 in char_table_2:
            sh = shared(ch1, ch2)
            count += sh
    return count


def quartet_distance(taxa: list[str], phylo1: PhyloTree, phylo2: PhyloTree):
    """If T1 and T2 are both unrooted binary trees on the same n taxa, then we now let q(T1,T2) denote the number of inferred quartets that are common to both trees. The quartet distance between T1 and T2, dq(T1,T2) is the number of quartets that are only inferred from one of the trees. More precisely, dq(T1,T2)=q(T1)+q(T2)-2q(T1,T2).

    Args:
        phylo1 (PhyloTree): _description_
        phylo2 (PhyloTree): _description_
    """
    # dq(T1, T1) = q(T1) + q(T2) - 2q(T1,T2)
    nleaves = len(taxa)
    qT1 = count_quartets(nleaves)
    qT2 = count_quartets(nleaves)
    dq = qT1 + qT2 - 2*count_shared_quartets2(taxa, phylo1, phylo2)

    return dq


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(taxa: list[str], newick1: str, newick2: str) -> int:
    phylo1 = PhyloTree(newick1)
    phylo2 = PhyloTree(newick2)
    dist = quartet_distance(taxa, phylo1, phylo2)
    return dist


def load_results(path: str) -> int:
    lines = readTextFile(path)
    count = int(lines[0])
    return count


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    taxa = lines[0].split(' ')
    newick1 = lines[1]
    newick2 = lines[2]

    count = solve(taxa, newick1, newick2)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(count, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_qrtd_0.txt'

    lines = readTextFile(path)
    taxa = lines[0].split(' ')
    newick1 = lines[1]
    newick2 = lines[2]

    # TODO: not implemented
    count = solve(taxa, newick1, newick2)

    # char_table_1 = set(phylo1.character_table_from_nontrivial_splits(taxa))
    # char_table_2 = set(phylo2.character_table_from_nontrivial_splits(taxa))

    # reduced = True
    # while reduced:
    #     reduced = False
    #     for ch1 in char_table_1:
    #         if ch1 in char_table_2:
    #             char_table_1.remove(ch1)
    #             char_table_2.remove(ch1)
    #             reduced = True
    #             break
    #         ch1_swap = swap_character_labels(ch1)
    #         if ch1_swap in char_table_2:
    #             char_table_1.remove(ch1)
    #             char_table_2.remove(ch1_swap)
    #             reduced = True
    #             break

    print(count)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(count), 'w')

    # correct = solve_and_check(path)
    # print(correct)
