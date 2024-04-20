from collections import Counter
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os


"""
https://rosalind.info/problems/asmq/

Given a collection of DNA strings representing contigs, we use the N statistic NXX (where XX ranges from 01 to 99) to represent the maximum positive integer L such that the total number of nucleotides of all contigs having length ≥L is at least XX% of the sum of contig lengths. The most commonly used such statistic is N50, although N75 is also worth mentioning.

    Given: A collection of at most 1000 DNA strings (whose combined length does not exceed 50 kbp).

    Return: N50 and N75 for this collection of strings.
"""

OutputT = tuple[int, int]

def Nstatistic(contigs: list[str], xx_target: int) -> int:
    assert 1 <= xx_target <= 99, "xx must be >=1 and <= 99"

    lengths = [len(contig) for contig in contigs]
    total = sum(lengths)
    counts = Counter(lengths)
    vals = sorted(list(set(lengths)), reverse=True)

    running_count = 0
    for l in vals:
        running_count += l*counts[l]
        xx_val = running_count/total*100
        if xx_val >= xx_target:
            return l
    return 0


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def solve(seqs: list[str]) -> OutputT:
    N50 = Nstatistic(seqs, 50)
    N75 = Nstatistic(seqs, 75)
    return (N50, N75)


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    N50, N75 = [int(v) for v in lines[0].split()]
    return (N50, N75)


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    seqs = [line for line in lines if len(line) and not line.isspace()]

    result = solve(seqs)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_asmq_1.txt'

    lines = readTextFile(path)
    seqs = [line for line in lines if len(line) and not line.isspace()]

    N50 = Nstatistic(seqs, 50)
    N75 = Nstatistic(seqs, 75)

    out = f"{N50} {N75}"
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
