
from BioInfoToolkit.Sequences.StringUtils import rotationally_equivalent

from BioInfoToolkit.Sequences.GenomeAssembly import DeBruijnMultiGraph, deBruijnMultiGraphFromReads
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
    https://rosalind.info/problems/grep/

Recall that a directed cycle is a cycle in a directed graph in which the head of one edge is equal to the tail of the following edge.

In a de Bruijn graph of k-mers, a circular string s is constructed from a directed cycle s1→s2→...→si→s1 is given by s1+s2[k]+...+si-k[k]+si-k+1[k]. That is, because the final k-1 symbols of s1 overlap with the first k-1 symbols of s2, we simply tack on the k-th symbol of s2 to s, then iterate the process.

For example, the circular string assembled from the cycle "AC" → "CT" → "TA" → "AC" is simply (ACT). Note that this string only has length three because the 2-mers "wrap around" in the string.

If every k-mer in a collection of reads occurs as an edge in a de Bruijn graph cycle the same number of times as it appears in the reads, then we say that the cycle is "complete."

    Given: A list Sk+1 of error-free DNA (k+1)-mers (k≤5) taken from the same strand of a circular chromosome (of length ≤50).

    Return: All circular strings assembled by complete cycles in the de Bruijn graph Bk of Sk+1. The strings may be given in any order, but each one should begin with the first (k+1)-mer provided in the input.
"""


def verify(result: list[str], solution: list[str]) -> bool:
    if len(result) != len(solution):
        return False
    sol_copy = solution.copy()
    for s1 in result:
        match = [s2 for s2 in sol_copy if rotationally_equivalent(s1, s2)]
        if not len(match):
            return False
        for m in match:
            sol_copy.remove(m)
        
    if len(sol_copy):
        return False
    return True


def solve(kmers: list[str]) -> list[str]:
    k1 = len(kmers[0])
    _g = deBruijnMultiGraphFromReads(kmers, k1)
    g = DeBruijnMultiGraph(_g, k1)
    # dot = g.draw_dot()
    # dot.view('g2')

    strings = [seq for seq in g.generateAllEulerianCycles(
        kmers[0][:-1], kmers[0][1:])]
    
    return strings

def load_results(path: str) -> list[str]:
    lines = readTextFile(path)
    strings = [line for line in lines if len(line) and not line.isspace()]
    return strings


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    string = solve(kmers)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(string, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_grep_1.txt'

    lines = readTextFile(path)
    kmers = [line for line in lines if len(line) and not line.isspace()]

    strings = solve(kmers)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')
    for s in strings:
        print(s)
        writeTextFile(result_path, s, 'a')

    correct = solve_and_check(path)
    print(correct)
