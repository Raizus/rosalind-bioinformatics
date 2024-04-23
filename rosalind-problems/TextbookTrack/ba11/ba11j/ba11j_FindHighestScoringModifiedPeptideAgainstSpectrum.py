from BioInfoToolkit.Spectrometry.Spectrometry import AminoacidMonoisotopicMassInteger, PSMGraph, findModifiedPeptide, parseModifiedPeptideToArray
import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile

"""
https://rosalind.info/problems/ba11j/

In addition to searching for mutated peptides, we will also need to search for post-translational modifications, which alter amino acids after a protein has been translated from RNA. In fact, most proteins are modified after translation, and hundreds of types of modifications have been discovered. For example, the enzymatic activity of many proteins is regulated by the addition or removal of a phosphate group at a specific amino acid. This process, called phosphorylation, is reversible; protein kinases add phosphate groups, whereas protein phosphatases remove them.

A modification of mass δ applied to an amino acid results in adding δ to the mass of this amino acid. For example, δ = 80 for phosphorylated amino acids (serine, threonine, and tyrosine), δ = 16 for the modification of proline into hydroxyproline, and δ = -1 for the modification of lysine into allysin. If δ is positive, then the resulting modified peptide has a peptide vector that differs from the original peptide vector Peptide' by inserting a block of δ zeroes before the i-th occurrence of 1 in Peptide'. In the more rare case that δ is negative, the modified peptide corresponds to deleting a block of |δ| zeroes from Peptide'.

We will use the term block indel to refer to the addition or removal of a block of consecutive zeroes from a binary vector. Thus, applying k modifications to an aminoacid string Peptide corresponds to applying k block indels to its peptide vector Peptide'. We define Variantsk(Peptide) as the set of all modified variants of Peptide with up to k modifications.

Given a peptide Peptide and a spectral vector Spectrum', our goal is to find a modified peptide from Variantsk(Peptide) with maximum score against Spectrum'.

Spectral Alignment Problem
    Given a peptide and a spectral vector, find a modified variant of this peptide that maximizes the peptide-spectrum score among all variants of the peptides with up to k modifications.

        Given: A peptide Peptide, a spectral vector Spectrum', and an integer k.

        Return: A peptide Peptide' related to Peptide by up to k modifications with maximal score against Spectrum' out of all possibilities.

Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty amino acids. Examples sometimes use imaginary amino acids X and Z having respective integer masses 4 and 5.
"""

OutputT = str

def verify(result: OutputT, solution: OutputT) -> bool:
    mass_array_res = parseModifiedPeptideToArray(result)
    mass_array_sol = parseModifiedPeptideToArray(solution)
    correct = mass_array_res == mass_array_sol
    return correct


def solve(peptide: str, spectral_vector: list[int], num_mods: int) -> OutputT:
    massDict = AminoacidMonoisotopicMassInteger

    # psm_graph = PSMGraph(peptide, spectral_vector, massDict)
    # d = psm_graph.draw_svg()
    # d.save_svg('example.svg')

    peptide_mod = findModifiedPeptide(peptide, spectral_vector, num_mods, massDict)
    return peptide_mod


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    mod_peptide = lines[0]

    return mod_peptide


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)

    peptide = lines[0]
    spectral_vector = [0] + [int(v) for v in lines[1].split()]
    num_mods = int(lines[2])

    result = solve(peptide, spectral_vector, num_mods)

    solution_path = solution_path_from_input_path(input_path)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba11j_1.txt'

    lines = readTextFile(path)

    peptide = lines[0]
    spectral_vector = [0] + [int(v) for v in lines[1].split()]
    num_mods = int(lines[2])

    mod_peptide = solve(peptide, spectral_vector, num_mods)

    out = str(mod_peptide)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
