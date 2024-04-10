
from collections import Counter, defaultdict
from itertools import combinations, product

def allPairwiseDifferences(L: list[int]):
    diffs: list[int] = []
    for a,b in product(L, repeat=2):
        diffs.append(a-b)

    return sorted(diffs)
    
def differenceMultiset(X: set[int]):
    diff: defaultdict[int, int] = defaultdict(int)
    for xi, xj in combinations(X, 2):
        diff[abs(xi-xj)] += 1
    return diff


def multisetDist(y: int, X: set[int]) -> Counter[int]:
    msDist = [abs(y-x) for x in X]
    return Counter(msDist)


def isSubset(a: Counter[int], b: Counter[int]) -> bool:
    """Returns True if a is a subset of b

    Args:
        a (Counter[int]): _description_
        b (Counter[int]): _description_

    Returns:
        bool: _description_
    """
    for key, value in a.items():
        if value > b.get(key, 0):
            return False
    return True


def multisetSubtract(a: Counter[int], b: Counter[int]) -> Counter[int]:
    # c = a - b
    c = a.copy()
    c.subtract(b)

    keys: list[int] = [key for key, val in c.items() if val <= 0]
    for key in keys:
        c.pop(key)

    return c


def multisetAdd(a: Counter[int], b: Counter[int]) -> Counter[int]:
    c = a + b
    return c


def setFromDiffMultiset(items: list[int]):
    """Solve the turnpike problem

    Args:
        items (list[int]): _description_

    Returns:
        _type_: _description_
    """
    items = sorted(items)
    width = items.pop()
    X = set([0, width])
    deltaX = Counter(items)

    possibleSets: list[set[int]] = []

    def place(deltaX: Counter[int], X: set[int]):
        # assume deltaX is sorted in ascending order
        if not len(deltaX):
            possibleSets.append(X.copy())
            return X
        y = max(deltaX.keys())
        distA = multisetDist(y, X)
        # if distA subset of deltaX
        if isSubset(distA, deltaX):
            X.add(y)
            deltaX = multisetSubtract(deltaX, distA)
            place(deltaX, X)
            X.remove(y)
            deltaX = multisetAdd(deltaX, distA)

        distB = multisetDist(width-y, X)
        # if distB subset of deltaX
        if isSubset(distB, deltaX):
            X.add(width-y)
            deltaX = multisetSubtract(deltaX, distB)
            place(deltaX, X)
            X.remove(width-y)
            deltaX = multisetAdd(deltaX, distB)

    place(deltaX, X)
    return possibleSets
    a = 0
