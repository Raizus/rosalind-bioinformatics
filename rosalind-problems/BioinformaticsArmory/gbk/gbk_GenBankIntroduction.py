from BioInfoToolkit.IO.IO import readTextFile
from Bio import Entrez
import re

# GenBankIntroduction
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_gbk.txt'
    path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_gbk.txt'

    lines = readTextFile(path)
    genus = lines[0]
    date1 = lines[1]
    date2 = lines[2]

    match1 = re.match(r'(\d{4})/(\d{1,2})/(\d{1,2})', date1)
    if match1:
        date1 = f"{int(match1[1]):04d}/{int(match1[2]):02d}/{int(match1[3]):02d}"
    match2 = re.match(r'(\d{4})/(\d{1,2})/(\d{1,2})', date2)
    if match2:
        date2 = f"{int(match2[1]):04d}/{int(match2[2]):02d}/{int(match2[3]):02d}"

    Entrez.email = "your_name@your_mail_server.com"
    term = f'"{genus}"[Organism]'

    handle = Entrez.esearch(db="nucleotide", term=term, datetype='pdat',
                            mindate=date1, maxdate=date2)
    record = Entrez.read(handle)
    handle.close()
    if record:
        print(record['Count']) # type: ignore

