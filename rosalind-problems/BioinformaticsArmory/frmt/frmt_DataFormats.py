from BioInfoToolkit.IO.IO import readTextFile
from Bio import Entrez
from Bio import SeqIO

# DataFormats
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_frmt.txt'
    path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_frmt.txt'

    lines = readTextFile(path)
    ids = lines[0].split()

    Entrez.email = "your_name@your_mail_server.com"
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype='fasta')
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    shortest = min(records, key = lambda x: len(x.seq))
    print(shortest.format("fasta"))
