import math
from BioInfoToolkit.IO.IO import readTextFile
from BioInfoToolkit.Sequences.StringUtils import suffixArray


def suffixArrayFindMatch(text: str, array: list[int], pattern: str) -> list[int]:
    i1, i2 = 0, len(array)

    matches: list[int] = []
    while i2 != i1:
        i_mid = math.floor((i2+i1)/2)
        pos = array[i_mid]
        suffix = text[pos:]
        if len(pattern) <= len(suffix) and pattern == suffix[:len(pattern)]:
            matches.append(pos)
        if suffix > pattern:
            i2 = i_mid
        else:
            i1 = i_mid
    return matches

# PatternMatchingWithSuffixArray
if __name__ == "__main__":
    path = 'rosalind-problems/TextbookTrack/test_data/test_ba9h.txt'
    # path = 'rosalind-problems/TextbookTrack/test_data/rosalind_ba9h.txt'

    lines = readTextFile(path)
    text = lines[0]
    patterns = lines[1:]

    text = text+'$'
    array = suffixArray(text)

    match1 = suffixArrayFindMatch(text, array, patterns[0])
    print(match1)


