import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation F_n = F_(n-1)+F_(n-2) and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months.

    Given: Positive integers n≤100 and m≤20.

    Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.
"""


def mortal_rabbits(n: int, m: int, f: int, C0: int = 1):
    """See https://arxiv.org/pdf/0710.2216.pdf
    Args:
        n (int): total months
        l (int): rabbits' life expectancy
        f (int): fertility age
        C0 (int): starting pop
        k (int): offspring count

    Returns:
        int: number of pairs of rabbits remaining after n months
    """
    Ci = C0
    pop: list[int] = [C0]
    for i in range(1, n):
        if i <= f: # rabbits haven't started reproducing, they'rejust reaching maturity
            Ci = Ci
        elif f + 1 <= i <= m + f - 2: # rabbits start reproducing but haven't started dying
            Ci = pop[i-1] + pop[i-f-1]
        else: # after rabbits start dying
            Ci = sum(pop[i-m-f+1:i-f])
        pop.append(Ci)
        print(f"Ci = {Ci}; i: {i}; f-1 = {f-1}; l + f - 2 = {m + f - 2}")
        print(f"i-l-f+1: {i-m-f+1}; i-f+1 = {i-f+1};")
        print(f"============================================")

    return pop



def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, m: int, f: int) -> int:
    pop = mortal_rabbits(n, m, f)
    return pop[-1]


def load_results(path: str) -> int:
    lines = readTextFile(path)
    pop = int(lines[0])
    return pop


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    (n, m) = [int(v) for v in lines[0].split()]
    f = 1

    result = solve(n, m, f)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_fibd_0.txt'

    lines = readTextFile(path)
    (n,m) = [int(v) for v in lines[0].split()]
    f = 1

    pop = solve(n,m,f)
    print(f"After {n} months there are {pop} rabbit pairs alive.")

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(pop), 'w')

