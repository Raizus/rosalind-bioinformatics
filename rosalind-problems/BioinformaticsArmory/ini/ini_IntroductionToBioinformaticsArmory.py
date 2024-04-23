from BioInfoToolkit.IO.IO import readTextFile
from BioInfoToolkit.Sequences.BioSequence import BioSequence
from BioInfoToolkit.Sequences.structures import DNA_NUCLEOTIDES


# IntroductionToBioinformaticsArmory
if __name__ == "__main__":
    path = 'rosalind-problems/BioinformaticsArmory/test_data/test_ini.txt'
    path = 'rosalind-problems/BioinformaticsArmory/test_data/rosalind_ini.txt'

    lines = readTextFile(path)
    sequence = lines[0]

    bioSeq = BioSequence(sequence)

    freq = bioSeq.nucleotide_frequency()
    vals = [freq[nuc] for nuc in DNA_NUCLEOTIDES]
    print(" ".join(str(val) for val in vals))