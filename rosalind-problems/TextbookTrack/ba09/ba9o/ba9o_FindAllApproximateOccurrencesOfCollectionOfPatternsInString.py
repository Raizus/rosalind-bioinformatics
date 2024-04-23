from BioInfoToolkit.IO.IO import readTextFile
from BioInfoToolkit.Sequences.BarrowsWheeler import BWTMultiplePatternMatching, barrowsWheelerTransform


def approximatePatternMatching(string: str, patterns: list[str], d: int):
    seeding = True

    bwt = barrowsWheelerTransform(string)

    if seeding:

        pass


# FindAllApproximateOccurrencesOfCollectionOfPatternsInString
if __name__ == "__main__":
    path = 'rosalind-problems/TextbookTrack/test_data/test_ba9o.txt'
    # path = 'rosalind-problems/TextbookTrack/test_data/rosalind_ba9o.txt'

    lines = readTextFile(path)
    string = lines[0]
    patterns = lines[1].split()
    d = int(lines[2])

