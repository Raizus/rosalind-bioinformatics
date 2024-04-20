from BioInfoToolkit.IO import getFastaDataFromProteinIds, readTextFile, writeTextFile
from BioInfoToolkit.Sequences.StringUtils import findMotif

import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows: [XY] means "either X or Y" and {X} means "any amino acid except X." For example, the N-glycosylation motif is written as N{P}[ST]{P}.

You can see the complete description and features of a particular protein by its access ID "uniprot_id" in the UniProt database, by inserting the ID number into

    http://www.uniprot.org/uniprot/uniprot_id

Alternatively, you can obtain a protein sequence in FASTA format by following

    http://www.uniprot.org/uniprot/uniprot_id.fasta

For example, the data for protein B5ZC00 can be found at http://www.uniprot.org/uniprot/B5ZC00.

    Given: At most 15 UniProt Protein Database access IDs.

    Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.
"""


def extract_id(header: str):
    return header.split('|')[1]


def verify(result: dict[str, list[int]], solution: dict[str, list[int]]) -> bool:
    for key, idxs in result.items():
        if key not in solution:
            return False
        if idxs != solution[key]:
            return False
    return True


def solve(fasta_dict: dict[str, str]) -> dict[str, list[int]]:
    pattern = "N[^P][ST][^P]"
    idxs_dict: dict[str, list[int]] = dict()
    for id, seq in fasta_dict.items():
        res = findMotif(pattern, seq)
        if len(res):
            indexes = [val+1 for val in res.keys()]
            idxs_dict[id] = indexes
    return idxs_dict


def load_results(path: str) -> dict[str, list[int]]:
    lines = readTextFile(path)

    idxs_dict: dict[str, list[int]] = dict()
    for i in range(len(lines)//2):
        id = lines[2*i]
        idxs = [int(v) for v in lines[2*i+1].split()]
        idxs_dict[id] = idxs

    return idxs_dict


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    id_label_dict: dict[str, str] = dict()
    for id in lines:
        _id = id.split("_")[0]
        id_label_dict[_id] = id

    print("Creating protein FASTA dictionary...")
    fasta_dict = getFastaDataFromProteinIds(list(id_label_dict.keys()))
    print("Dictionary created.")

    idxs_dict = solve(fasta_dict)
    idxs_dict = {id_label_dict[key]: val for key, val in idxs_dict.items()}

    solution = load_results(solution_path)

    correct = verify(idxs_dict, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_mprt_1.txt'

    lines = readTextFile(path)
    id_label_dict: dict[str, str] = dict()
    for id in lines:
        _id = id.split("_")[0]
        id_label_dict[_id] = id

    print("Creating protein FASTA dictionary...")
    fasta_dict = getFastaDataFromProteinIds(list(id_label_dict.keys()))
    print("Dictionary created.")

    idxs_dict = solve(fasta_dict)
    idxs_dict = {id_label_dict[key]: val for key, val in idxs_dict.items()}

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, None, 'w')

    for id, idxs in idxs_dict.items():
        writeTextFile(result_path, id, 'a')
        out = ' '.join([str(v) for v in idxs])
        writeTextFile(result_path, out, 'a')

        print(id)
        print(out)
