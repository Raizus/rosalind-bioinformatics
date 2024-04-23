from BioInfoToolkit.IO.IO import readTextFile, writeTextFile
from Bio import Entrez
from Bio import SeqIO

# FASTQ format introduction
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_tfsq.txt'
    path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_tfsq.txt'

    lines = readTextFile(path)

    records = list(SeqIO.parse(path, format='fastq'))

    path_res = "rosalind-problems/BioinformaticsArmory/results_files/tfsq_result.txt"
    writeTextFile(path_res, None, 'w')
    for record in records:
        print(record.format('fasta'))
        writeTextFile(path_res, record.format('fasta'), 'a')

