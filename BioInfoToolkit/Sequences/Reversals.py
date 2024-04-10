# https://archive.is/0322E

def performReversal(seq: list[int|str], start_idx: int, end_idx: int):
    prefix = seq[:start_idx]
    reversed_subseq = seq[start_idx:end_idx][::-1]
    suffix = seq[end_idx:]
    return prefix + reversed_subseq + suffix


def findBreakpoints(sequence: list[int | str], target: list[int | str]) -> list[int]:
    # a = 1, 2…4, 3…5, 6, 7, 8, 9, 10
    # b = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

    n = len(sequence)
    breakpoints: list[int] = []
    for idx in range(n-1):
        curr_el = sequence[idx]
        next_el = sequence[idx+1]

        target_curr_idx = target.index(curr_el)
        target_next_idx = target.index(next_el)

        if abs(target_curr_idx - target_next_idx) != 1:
            breakpoints.append(idx+1)
    return breakpoints


def findMinimumBreakpointReversals(sequences: list[list[int | str]], target: list[int | str]):
    reversals: list[list[int | str]] = []

    for seq in sequences:
        breakpoints = findBreakpoints(seq, target)
        for start_idx in range(len(breakpoints)-1):
            for end_idx in range(1, len(breakpoints)):
                reversal = performReversal(
                    seq, breakpoints[start_idx], breakpoints[end_idx])
                reversals.append(reversal)

    min_bp = len(target)
    minimum_reversals = []
    for reversal in reversals:
        num_breakpoints = len(findBreakpoints(reversal, target))
        if num_breakpoints < min_bp:
            min_bp = num_breakpoints
            minimum_reversals = [reversal]
        elif num_breakpoints == min_bp:
            minimum_reversals.append(reversal)
    return minimum_reversals


def findMinimumBreakpointReversalsWithHistories(sequences: list[tuple[list[int | str], list[tuple[int, int]]]], target: list[int | str]):
    reversalsHistory: list[tuple[list[int | str], list[tuple[int, int]]]] = []
    for seq in sequences:
        breakpoints = findBreakpoints(seq[0], target)
        for start_idx in range(len(breakpoints)-1):
            for end_idx in range(1, len(breakpoints)):
                reversal = performReversal(
                    seq[0], breakpoints[start_idx], breakpoints[end_idx])

                reversalHist = (reversal, seq[1] + [(breakpoints[start_idx]-1, breakpoints[end_idx]-1)])
                reversalsHistory.append(reversalHist)
    
    min_bp = len(target)
    minimum_reversals = []

    for reversal in reversalsHistory:
        num_breakpoints = len(findBreakpoints(reversal[0], target))
        if num_breakpoints < min_bp:
            min_bp = num_breakpoints
            minimum_reversals = [reversal]
        elif num_breakpoints == min_bp:
            minimum_reversals.append(reversal)
    return minimum_reversals

def reversalDistance(sequence: list[int], target: list[int]) -> int:
    n1 = len(sequence)
    n2 = len(target)
    assert n1 == n2

    padded_sequence = ['-'] + sequence + ['+']
    padded_target = ['-'] + target + ['+']

    reversals = 0
    current_sequences = [padded_sequence]

    while padded_target not in current_sequences:
        current_sequences = findMinimumBreakpointReversals(
            current_sequences, padded_target)
        reversals += 1
    return reversals


def getReversalDistanceWithHistories(sequence: list[int], target: list[int]):
    padded_sequence = ['-'] + sequence + ['+']
    padded_target = ['-'] + target + ['+']

    reversals = 0
    current_sequences: list[tuple[list[int | str], list[tuple[int, int]]]] = [
        (padded_sequence, [])]
    while padded_target not in [current_sequence[0] for current_sequence in current_sequences]:
        current_sequences = findMinimumBreakpointReversalsWithHistories(
            current_sequences, padded_target)
        reversals += 1
    return reversals, current_sequences