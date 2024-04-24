from BioInfoToolkit.IO.IO import readTextFile, writeTextFile


if __name__ == "__main__":
    path = 'rosalind-problems/TextbookTrack/test_data/test_ba6a.txt'
    # path = 'rosalind-problems/TextbookTrack/test_data/rosalind_ba6a.txt'

    lines = readTextFile(path)
    line = lines[0]
    values = [v for v in line.rstrip(')').lstrip('(').split()]
    print(values)

