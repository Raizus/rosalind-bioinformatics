
from time import sleep
import warnings
import requests
import os
import re

def readTextFile(filePath: str) -> list[str]:
    with open(filePath, 'r') as f:
        return [l.strip() for l in f.readlines()]


def writeTextFile(filePath: str, text: str | None, mode="w"):
    with open(filePath, mode) as f:
        if text is None:
            f.write('')
        else:
            f.write(text + '\n')


def FASTA_lines_to_dict(lines: list[str]):
    lines = [l.strip() for l in lines]
    FASTADict: dict[str, str] = {}
    FASTALabel = ""
    for line in lines:
        if '>' in line:
            FASTALabel = line[1:]
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line

    return FASTADict


def read_FASTA(filePath: str):
    with open(filePath, 'r') as f:
        FASTAData = [l.strip() for l in f.readlines()]

    return FASTA_lines_to_dict(FASTAData)


def getFastaDataFromProteinIds(cIDs: str | list[str]):
    baseUrl = "http://www.uniprot.org/uniprot/"
    sleepT = 0.2

    if isinstance(cIDs, str):
        cIDs = [cIDs]

    FASTADict: dict[str, str] = {}
    for cID in cIDs:
        currentUrl = baseUrl+cID+".fasta"
        response = requests.post(currentUrl)
        if not response.ok:
            warnings.warn(f"Protein with id {cID} not found. Skipping...")
            sleep(sleepT)
            continue

        FASTAData = response.text.splitlines()
        FASTALabel = ""
        for line in FASTAData:
            if '>' in line:
                FASTALabel = line[1:]
                FASTADict[FASTALabel] = ""
            else:
                FASTADict[FASTALabel] += line.strip()
        sleep(sleepT)

    return FASTADict


def get_input_file_suffix(path: str) -> str:
    tail, head = os.path.split(path)
    filename, extension = os.path.splitext(head)

    match = re.match(r'rosalind_(\w+(?:_(?:\d+)))', filename)
    if match is None:
        raise Exception("File name doesn't match template. File name must have the format 'rosalind_{suffix}.txt'.")
    suffix = match[1]

    return suffix


def solution_path_from_input_path(input_path: str) -> str:
    suffix = get_input_file_suffix(input_path)
    tail, head = os.path.split(input_path)
    solution_path = f'{tail}/solution_{suffix}.txt'
    return solution_path


def result_path_from_input_path(input_path: str) -> str:
    suffix = get_input_file_suffix(input_path)
    tail, head = os.path.split(input_path)
    result_path = f'{tail}/result_{suffix}.txt'
    return result_path
