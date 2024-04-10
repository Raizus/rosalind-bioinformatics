
from typing import Iterable


def generateStringCycles(string: str):
    for i in range(len(string)):
        yield string[i:] + string[:i]


def barrowsWheelerTransform(string: str) -> str:
    sortedStrings = sorted(s for s in generateStringCycles(string))
    bwt = ''.join([s[-1] for s in sortedStrings])
    return bwt


def index_strings_by_occurrence(array: Iterable[str]) -> list[tuple[str, int]]:
    """Given a list of strings retuns a list of tuples where each string is associated with an index
    corresponding to the k-th occurence of that string in the list (0-indexed)

    Args:
        array (list[str]): _description_

    Returns:
        list[tuple[str, int]]: _description_
    """
    counts: dict[str, int] = dict()
    result: list[tuple[str, int]] = []
    for char in array:
        count = counts.get(char, 0)
        result.append((char, count))
        counts[char] = count + 1
    return result


def stringFromBWT(bwt: str) -> str:
    """Memory efficient string reconstrution from Barrows-Wheeler Transform

    Args:
        bwt (str): _description_

    Returns:
        str: reconstructed string
    """
    last_column = [c for c in bwt]
    first_column = sorted(last_column)    
    last_column = index_strings_by_occurrence(last_column)
    first_column = index_strings_by_occurrence(first_column)

    temp_dict = {a: b for a, b in zip(first_column, last_column)}

    current = first_column[0]
    string = current[0]
    for i in range(len(bwt)-1):
        next_symbol = temp_dict[current]
        string = next_symbol[0] + string
        current = next_symbol

    return string


def lastToFirstArray(bwt: str) -> list[int]:
    """Given a symbol at position i in LastColumn (the bwt), what is its position in FirstColumn (the fisrt symbol in the list of lexicographically sorted cyclic strings)

    Args:
        bwt (str): _description_

    Returns:
        list[int]: last to first array
    """
    first_to_last = [i[0] for i in sorted(enumerate(bwt), key=lambda x: x[1])]
    last_to_first = [0 for _ in bwt]
    for i, v in enumerate(first_to_last):
        last_to_first[v] = i
    return last_to_first


def lastToFirst(bwt: str, i: int) -> int:
    """Given a symbol at position i in LastColumn (the bwt), what is its position in FirstColumn (the first symbol in the list of lexicographically sorted cyclic strings)

    Args:
        bwt (str): _description_
        i (int): position of symbol i in the bwt

    Returns:
        int: position k in the last first column, of symbol i in first column. 
            k: bwt[i] == first_column[k]
    """
    first_to_last = [j[0] for j in sorted(enumerate(bwt), key=lambda x: x[1])]
    return first_to_last.index(i)


def suffixArrayFromBWT(bwt: str, last_to_first: list[int] | None = None):
    l = len(bwt)

    if last_to_first is None:
        last_to_first = lastToFirstArray(bwt)

    suffix_array: list[int] = [0 for _ in bwt]
    i = 0
    si = l - 1
    suffix_array[i] = si
    while si > 0:
        i = last_to_first[i]
        si -= 1
        suffix_array[i] = si

    return suffix_array


def buildCheckpointCountArraysBWT(bwt: str, c: int):
    l = len(bwt)

    checkpoint_array: list[dict[str, int]] = []
    counts: dict[str, int] = dict()
    checkpoint_array.append(dict())

    for i, char in enumerate(bwt, 1):
        counts.setdefault(char, 0)
        counts[char] += 1
        if i % c == 0:
            checkpoint_array.append(counts)
            counts = counts.copy()
    return checkpoint_array


class CheckpointArrayBWT:
    array: list[dict[str, int]]
    c: int

    def __init__(self, bwt: str, c: int) -> None:
        """Creates a checkpoint array for BWTMatching

        Args:
            bwt (str): BWT string
            c (int): creates checkpoints when i % c == 0 (i is the index of the bwt string)
        """
        l = len(bwt)

        checkpoint_array: list[dict[str, int]] = []
        counts: dict[str, int] = dict()
        checkpoint_array.append(dict())

        for i, char in enumerate(bwt, 1):
            counts.setdefault(char, 0)
            counts[char] += 1
            if i % c == 0:
                checkpoint_array.append(counts)
                counts = counts.copy()

        self.array = checkpoint_array
        self.c = c

    def get_count(self, bwt: str, symbol: str, i: int):
        c = self.c
        d = min(round(i / c), len(bwt)//c)
        count = self.array[d].get(symbol, 0)
        if d <= i / c:
            for i in range(d*c, i+1):
                if bwt[i] == symbol:
                    count += 1
        else:
            for i in range(d*c-1, i, -1):
                if bwt[i] == symbol:
                    count -= 1
        return count


def getFirstOccurences(first_column: list[str]) -> dict[str, int]:
    first_occurence: dict[str, int] = dict()
    for i, symbol in enumerate(first_column):
        first_occurence.setdefault(symbol, i)

    return first_occurence


def BWTMatchingWithCheckpoints(bwt: str, pattern: str, c: int = 5,
                               checkpoint_array: CheckpointArrayBWT | None = None,
                               first_occurence: dict[str, int] | None = None):
    """Performs Borrows-Wheeler Transform (BWT) matching, using checkpoint arrays. When performing multiple matches, providing a pre-computed checkpoint array and first ocurrences dictionary will speed up performance significantly

    Args:
        bwt (str): the BWT string
        pattern (str): the pattern to match
        c (int, optional): creates checkpoints when i % c == 0 (i is the index of the bwt string), when checkpoint_array is not provided. Defaults to 5.
        checkpoint_array (CheckpointArrayBWT | None, optional): Checkpoint array. Defaults to None.
        first_occurence (dict[str, int] | None, optional): Dictionary of the first ocurrence of each symbol in the first column of the BWT. Defaults to None.

    Returns:
        _type_: _description_
    """
    last_column = [char for char in bwt]
    first_column = sorted(last_column)

    if checkpoint_array is None:
        checkpoint_array = CheckpointArrayBWT(bwt, c)

    # first occurence in first column
    if first_occurence is None:
        first_occurence = getFirstOccurences(first_column)

    top = 0
    bottom = len(last_column)-1
    i = len(pattern)-1
    while top <= bottom:
        if i >= 0:
            symbol = pattern[i]
            i -= 1
            if symbol in last_column[top:bottom+1]:
                # topIdx = idx of first ocurrence of symbol from top to bottom in last_column
                # bottomIdx = idx of the last ocurrence of symbol
                # top = lastToFirst(topIndex)
                # bottom = lastToFirst(bottomIndex)
                top = first_occurence[symbol] + checkpoint_array.get_count(bwt, symbol, top-1)
                bottom = first_occurence[symbol] + checkpoint_array.get_count(bwt, symbol, bottom) - 1
            else:
                return None
        else:
            return top, bottom


def BWTMultiplePatternMatching(bwt: str, patterns: list[str], c: int):
    """_summary_

    Args:
        bwt (str): the Barrows-Wheeler Transform
        patterns (list[str]): the patterns to match
        c (int): creates checkpoints when i % c == 0 (i is the index of the bwt string), when checkpoint_array is not provided.

    Returns:
        matches list[int]: array with the match index position in the original string
    """
    checkpoint_array = CheckpointArrayBWT(bwt, c)
    first_occurence = getFirstOccurences(sorted(char for char in bwt))

    suffix_array = suffixArrayFromBWT(bwt)

    matches: list[int] = []
    for pattern in patterns:
        match = BWTMatchingWithCheckpoints(bwt, pattern, c, checkpoint_array, first_occurence)
        if match is not None:
            top, bottom = match
            matches.extend(suffix_array[i] for i in range(top, bottom+1))

    return matches


def BWTInexactMatchingWithCheckpoints(bwt: str, pattern: str, d: int, c: int = 5,
                               checkpoint_array: CheckpointArrayBWT | None = None,
                               first_occurence: dict[str, int] | None = None):
    """Performs Borrows-Wheeler Transform (BWT) matching, using checkpoint arrays. When performing multiple matches, providing a pre-computed checkpoint array and first ocurrences dictionary will speed up performance significantly

    Args:
        bwt (str): the BWT string
        pattern (str): the pattern to match
        c (int, optional): creates checkpoints when i % c == 0 (i is the index of the bwt string), when checkpoint_array is not provided. Defaults to 5.
        checkpoint_array (CheckpointArrayBWT | None, optional): Checkpoint array. Defaults to None.
        first_occurence (dict[str, int] | None, optional): Dictionary of the first ocurrence of each symbol in the first column of the BWT. Defaults to None.

    Returns:
        _type_: _description_
    """
    last_column = [char for char in bwt]
    first_column = sorted(last_column)

    mismatch_counts = [0 for _ in bwt]

    if checkpoint_array is None:
        checkpoint_array = CheckpointArrayBWT(bwt, c)

    # first occurence in first column
    if first_occurence is None:
        first_occurence = getFirstOccurences(first_column)


    for i in range(len(pattern)-1, -1, -1):
        symbol = pattern[i]

    top = 0
    bottom = len(last_column)-1
    i = len(pattern)-1
    while top <= bottom:
        if i >= 0:
            symbol = pattern[i]
            i -= 1
            if symbol in last_column[top:bottom+1]:
                # topIdx = idx of first ocurrence of symbol from top to bottom in last_column
                # bottomIdx = idx of the last ocurrence of symbol
                # top = lastToFirst(topIndex)
                # bottom = lastToFirst(bottomIndex)
                top = first_occurence[symbol] + \
                    checkpoint_array.get_count(bwt, symbol, top-1)
                bottom = first_occurence[symbol] + \
                    checkpoint_array.get_count(bwt, symbol, bottom) - 1
            else:
                return None
        else:
            return top, bottom

