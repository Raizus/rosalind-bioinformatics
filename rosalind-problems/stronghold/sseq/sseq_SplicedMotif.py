
from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, result_path_from_input_path, writeTextFile
import os

"""
A subsequence of a string is a collection of symbols contained in order (though not necessarily contiguously) in the string (e.g., ACG is a subsequence of TATGCTAAGATC). The indices of a subsequence are the positions in the string at which the symbols of the subsequence appear; thus, the indices of ACG in TATGCTAAGATC can be represented by (2, 5, 9).

As a substring can have multiple locations, a subsequence can have multiple collections of indices, and the same index can be reused in more than one appearance of the subsequence; for example, ACG is a subsequence of AACCGGTT in 8 different ways.

    Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.

    Return: One collection of indices of s in which the symbols of t appear as a subsequence of s. If multiple solutions exist, you may return any one.
"""


def buildPattern(seq: str) -> str:
    pattern = ".*".join(f"({s})" for s in seq)
    return pattern


def find_subsequence_motif(seq: str, motif: str) -> list[int]:
    i1 = 0
    s1 = motif[i1]
    idxs: list[int] = []
    for i2, s2 in enumerate(seq):
        if s1 != s2:
            continue
        idxs.append(i2)
        i1 += 1
        if i1 >= len(motif):
            break
        s1 = motif[i1]
    if len(idxs) != len(motif):
        raise Exception('Motif is not a subsequence of seq.')
    return idxs


def verify(result: list[int], seq: str, motif: str) -> bool:
    subseq = ''.join(seq[r-1] for r in result)
    strict_increasing = all(v1<v2 for v1, v2 in zip(result, result[1:]))
    correct = result[-1] <= len(seq) and subseq == motif and strict_increasing
    return correct


def solve(seq: str, motif: str) -> list[int]:
    idxs = find_subsequence_motif(seq, motif)
    idxs = [i+1 for i in idxs]
    # pattern = buildPattern(motif)
    # pattern = re.compile(f"(?=({pattern}))")

    # idxs: list[int] = []
    # for m in pattern.finditer(seq):
    #     idxs = [m.start(i)+1 for i in range(2, len(m.regs))]
    #     return idxs
    return idxs

def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    idxs = [int(v) for v in lines[0].split()]
    return idxs


def solve_and_check(input_path: str) -> bool:
    fasta_dict = read_FASTA(input_path)
    seqs = sorted(fasta_dict.values(), key=lambda s: len(s), reverse=True)
    seq = seqs[0]
    motif = seqs[1]
    result = solve(seq, motif)

    # solution_path = solution_path_from_input_path(input_path)
    # solution = load_results(solution_path)

    correct = verify(result, seq, motif)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_sseq_1.txt'

    fasta_dict = read_FASTA(path)
    seqs = sorted(fasta_dict.values(), key=lambda s: len(s), reverse=True)
    seq = seqs[0]
    motif = seqs[1]

    idxs = solve(seq, motif)

    out = ' '.join(str(idx) for idx in idxs)
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    # correct = solve_and_check(path)
    # print(correct)