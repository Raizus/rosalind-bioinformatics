from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

from BioInfoToolkit.Sequences.SequenceUtils import gc_content

"""
The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.

    Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

    Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.
"""


def verify(result: tuple[str, float], solution: tuple[str, float]) -> bool:
    error = abs(result[1] - solution[1])
    correct = result[0] == solution[0] and error <= 0.001

    return correct


def solve(fasta_dict: dict[str, str]) -> tuple[str, float]:
    gc_content_dict = {label: gc_content(seq) for label, seq in fasta_dict.items()}
    max_gc_key = max(gc_content_dict, key=lambda label: gc_content_dict[label])
    result = (max_gc_key, gc_content_dict[max_gc_key]*100)
    return result

def load_results(path: str) -> tuple[str, float]:
    lines = readTextFile(path)
    label = lines[0]
    gc_value = float(lines[1])
    return (label, gc_value)


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    fasta_dict = read_FASTA(input_path)
    result = solve(fasta_dict)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_gc_1.txt'

    fasta_dict = read_FASTA(path)
    result = solve(fasta_dict)
    print(f"{result[0]}\n{result[1]}")

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, result[0], 'w')
    writeTextFile(result_path, str(result[1]), 'a')
