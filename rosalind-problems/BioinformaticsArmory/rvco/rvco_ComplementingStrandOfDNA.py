from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, writeTextFile
from Bio import Entrez
from Bio import SeqIO

from BioInfoToolkit.Sequences.BioSequence import BioSequence

# ComplementingStrandOfDNA
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_rvco.txt'
    path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_rvco.txt'

    fasta_dict = read_FASTA(path)

    matches_count = sum([BioSequence(seq).reverseComplement() == seq for seq in fasta_dict.values()])
    print(matches_count)