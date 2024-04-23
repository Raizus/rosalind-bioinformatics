from BioInfoToolkit.IO.IO import read_FASTA, readTextFile, writeTextFile
from Bio import Entrez
from Bio import SeqIO

from BioInfoToolkit.Sequences.BioSequence import BioSequence

# FindingGenesWithORFs
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_orfr.txt'
    # path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_orfr.txt'

    lines = readTextFile(path)
    seq = lines[0]

    bioSeq = BioSequence(seq, 'DNA')

    all_prots = bioSeq.all_proteins_from_rfs()
    max_length = max([len(prot) for prot in all_prots])
    prots = [prot for prot in all_prots if len(prot) == max_length]
    for prot in prots:
        print(prot)
