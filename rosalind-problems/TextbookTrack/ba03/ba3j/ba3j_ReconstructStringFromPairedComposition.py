from BioInfoToolkit.Sequences.GenomeAssembly import PairedDeBruijnMultiGraph, deBruijnMultiGraphFromPairedKmers
from BioInfoToolkit.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/ba3j/

Since increasing read length presents a difficult experimental problem, biologists have suggested an indirect way of increasing read length by generating read-pairs, which are pairs of reads separated by a fixed distance d in the genome.

You can think about a read-pair as a long "gapped" read of length k + d + k whose first and last k-mers are known but whose middle segment of length d is unknown. Nevertheless, read-pairs contain more information than k-mers alone, and so we should be able to use them to improve our assemblies. If only you could infer the nucleotides in the middle segment of a read-pair, you would immediately increase the read length from k to 2 · k + d.

Given a string Text, a (k,d)-mer is a pair of k-mers in Text separated by distance d. We use the notation (Pattern1|Pattern2) to refer to a a (k,d)-mer whose k-mers are Pattern1 and Pattern2. The (k,d)-mer composition of Text, denoted PairedCompositionk,d(Text), is the collection of all (k,d)- mers in Text (including repeated (k,d)-mers).

String Reconstruction from Read-Pairs Problem
    Reconstruct a string from its paired composition.

        Given: Integers k and d followed by a collection of paired k-mers PairedReads.

        Return: A string Text with (k, d)-mer composition equal to PairedReads. (If multiple answers exist, you may return any one.)
"""

OutputT = str


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution

    return correct


def solve(pairs: list[tuple[str,str]], k: int, d: int) -> OutputT:
    _g = deBruijnMultiGraphFromPairedKmers(pairs)
    g = PairedDeBruijnMultiGraph(_g, k, d+k)
    # dot = g.draw_dot()
    # dot.view()

    string = g.reconstructStringFromEulerianPath()
    return string


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    sequence = lines[0]
    return sequence


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(input_path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, d = [int(v) for v in lines[0].split()]
    pairs = [tuple(a.split('|')[:2]) for a in lines[1:]]
    pairs = [(pair[0], pair[1]) for pair in pairs]

    result = solve(pairs, k, d)

    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba3j_1.txt'

    lines = readTextFile(path)
    lines = [line for line in lines if len(line) and not line.isspace()]
    k, d = [int(v) for v in lines[0].split()]
    pairs = [tuple(a.split('|')[:2]) for a in lines[1:]]
    pairs = [(pair[0], pair[1]) for pair in pairs]

    string = solve(pairs, k, d)

    out = string
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)

