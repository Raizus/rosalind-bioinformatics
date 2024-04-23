from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, writeTextFile
from Bio import Entrez
from Bio import SeqIO

from BioInfoToolkit.Sequences.BioSequence import BioSequence

# NewMotifDiscovery
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_meme.txt'
    # path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_meme.txt'

