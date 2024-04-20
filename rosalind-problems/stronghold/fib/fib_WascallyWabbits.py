import os
from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile


"""
https://rosalind.info/problems/fib/

A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence (π,-√2,0,π) and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation an to represent the n-th term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fn=Fn-1+Fn-2 (with F1=F2=1 to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.

    Given: Positive integers n≤40 and k≤5.

    Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
"""


def rabbits(generation: int, offspring_count: int) -> int:
    # r[n] = k*r[n-1] = k*r[n-1]
    parent, child = 1, 1
    for _ in range(generation - 1):
        child, parent = parent, parent + child * offspring_count
    return child


def verify(result: int, solution: int) -> bool:
    correct = result == solution
    return correct


def solve(n: int, m: int) -> int:
    pop = rabbits(n, m)
    return pop


def load_results(path: str) -> int:
    lines = readTextFile(path)
    pop = int(lines[0])
    return pop


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    (n, m) = [int(v) for v in lines[0].split()]

    result = solve(n, m)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_fib_1.txt'

    lines = readTextFile(path)
    (n, m) = [int(v) for v in lines[0].split()]

    pop = solve(n, m)
    print(f"After {n} months there are {pop} rabbit pairs alive.")

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(pop), 'w')

    correct = solve_and_check(path)
    print(correct)
