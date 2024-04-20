from BioInfoToolkit.Sequences.GenomeAssembly import overlap_graph_assembly, superstring_from_overlap_graph
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
For a collection of strings, a larger string containing every one of the smaller strings as a substring is called a superstring.

By the assumption of parsimony, a shortest possible superstring over a collection of reads serves as a candidate chromosome.

    Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format (which represent reads deriving from the same strand of a single linear chromosome).
    The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

    Return: A shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).
"""


def verify(result: str, solution: str) -> bool:
    return result == solution


def solve(fasta_dict: dict[str,str]) -> str:
    graph = overlap_graph_assembly(fasta_dict)
    superstring = superstring_from_overlap_graph(fasta_dict, graph)
    return superstring


def load_results(path: str) -> str:
    lines = readTextFile(path)
    superstring = lines[0]
    return superstring


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(path)
    result = solve(fasta_dict)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_long_1.txt'

    fasta_dict = read_FASTA(path)
    graph = overlap_graph_assembly(fasta_dict)
    superstring = superstring_from_overlap_graph(fasta_dict, graph)

    print(superstring)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, superstring, 'w')

    correct = solve_and_check(path)
    print(correct)