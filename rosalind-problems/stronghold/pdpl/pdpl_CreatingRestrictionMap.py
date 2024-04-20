from collections import Counter
from BioInfoToolkit.Sets import differenceMultiset, setFromDiffMultiset

from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
https://rosalind.info/problems/pdpl/

For a set X containing numbers, the difference multiset of X is the multiset ΔX defined as the collection of all positive differences between elements of X. As a quick example, if X={2,4,7}, then we will have that ΔX={2,3,5}.

If X contains n elements, then ΔX will contain one element for each pair of elements from X, so that ΔX contains (n2) elements (see combination statistic). You may note the similarity between the difference multiset and the Minkowski difference X⊖X, which contains the elements of ΔX and their negatives. For the above set X, X⊖X is {-5,-3,-2,2,3,5}.

In practical terms, we can easily obtain a multiset L corresponding to the distances between restriction sites on a chromosome. If we can find a set X whose difference multiset ΔX is equal to L, then X will represent possible locations of these restriction sites. For an example, consult Figure 1.

    Given: A multiset L containing (n2) positive integers for some positive integer n.

    Return: A set X containing n nonnegative integers such that ΔX=L.
"""

# https://www.cs.ucf.edu/courses/cap5510/fall2009/res.map/Restrict.Mapping.pdf

def verify(result: list[int], solution: list[int], ms: list[int]) -> bool:
    diff = differenceMultiset(set(result))
    correct = Counter(ms) == Counter(diff)
    return correct


def solve(ms: list[int]) -> list[int]:
    possible_sets = setFromDiffMultiset(ms)
    return sorted(possible_sets[0])


def load_results(path: str) -> list[int]:
    lines = readTextFile(path)
    set = [int(v) for v in lines[0].split()]
    return set


def solve_and_check(input_path: str) -> bool:
    lines = readTextFile(input_path)
    values = [int(val) for val in lines[0].split()]

    possible_set = solve(values)

    solution_path = solution_path_from_input_path(input_path)
    solution = load_results(solution_path)

    correct = verify(possible_set, solution, values)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_pdpl_1.txt'

    lines = readTextFile(path)
    values = [int(val) for val in lines[0].split()]

    possible_sets = setFromDiffMultiset(values)

    out = ' '.join(str(val) for val in possible_sets[0])
    print(out)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    correct = solve_and_check(path)
    print(correct)
