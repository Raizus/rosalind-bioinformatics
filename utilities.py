
NUCLEOTIDE_COLORMAP = {
    'A': "\033[92m",
    'C': "\033[94m",
    'G': "\033[93m",
    'T': "\033[91m",
    'U': "\033[91m",
    'reset': "\033[0;0m",
}


def colored(seq: str) -> str:
    tempStr = ""
    for nuc in seq:
        if nuc in NUCLEOTIDE_COLORMAP:
            tempStr += NUCLEOTIDE_COLORMAP[nuc] + nuc
        else:
            tempStr += NUCLEOTIDE_COLORMAP['reset'] + nuc

    return tempStr + "\033[0;0m"
