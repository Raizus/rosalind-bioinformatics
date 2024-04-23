from BioInfoToolkit.IO.IO import readTextFile, writeTextFile
from Bio import Entrez
from Bio import SeqIO

from BioInfoToolkit.Sequences.BioSequence import BioSequence

# ProteinTranslation
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_ptra.txt'
    path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_ptra.txt'

    lines = readTextFile(path)
    seq = lines[0]
    prot = lines[1]

    bioSeq = BioSequence(seq, 'DNA')
    transl = ''.join(bioSeq.translate_seq())
    trans = transl.maketrans({'_': 'Q'})
    transl = transl.translate(trans)
    print(transl)
    idx = transl.find(prot)
    if idx != -1:
        print(idx+1)


