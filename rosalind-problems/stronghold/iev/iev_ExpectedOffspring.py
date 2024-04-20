from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
import os

"""
For a random variable X taking integer values between 1 and n, the expected value of X is E(X)=∑nk=1k×Pr(X=k). The expected value offers us a way of taking the long-term average of a random variable over a large number of trials.

As a motivating example, let X be the number on a six-sided die. Over a large number of rolls, we should expect to obtain an average of 3.5 on the die (even though it's not possible to roll a 3.5). The formula for expected value confirms that E(X)=∑6k=1k×Pr(X=k)=3.5.

More generally, a random variable for which every one of a number of equally spaced outcomes has the same probability is called a uniform random variable (in the die example, this "equal spacing" is equal to 1). We can generalize our die example to find that if X is a uniform random variable with minimum possible value a and maximum possible value b, then E(X)=a+b2. You may also wish to verify that for the dice example, if Y is the random variable associated with the outcome of a second die roll, then E(X+Y)=7.

    Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the following genotypes:
        1. AA-AA
        2. AA-Aa
        3. AA-aa
        4. Aa-Aa
        5. Aa-aa
        6. aa-aa

    Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.
"""

def expected_offspring_dominant_phenotype(genotype_counts: list[int]) -> float:
    """Returns the expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.

    Args:
        genotype_counts (list[int]): 
            1. AA-AA
            2. AA-Aa
            3. AA-aa
            4. Aa-Aa
            5. Aa-aa
            6. aa-aa


    Returns:
        float: expected population
    """
    k = 2
    expected = genotype_counts[0]*k + genotype_counts[1]*k + genotype_counts[2]*k \
        + genotype_counts[3]*k*3/4 + genotype_counts[4] * \
        k*1/2 + genotype_counts[5]*k*0
    return expected


def verify(result: float, solution: float) -> bool:
    error = abs(result - solution)
    correct = error <= 0.001

    return correct


def solve(genotype_pairs_count: list[int]) -> float:
    result = expected_offspring_dominant_phenotype(genotype_pairs_count)
    return result


def load_results(path: str) -> float:
    lines = readTextFile(path)
    value = float(lines[0])
    return value


def solve_and_check(input_path: str) -> bool:
    solution_path = solution_path_from_input_path(input_path)

    lines = readTextFile(path)
    genotype_pairs_count = [int(v) for v in lines[0].split()]
    result = solve(genotype_pairs_count)
    solution = load_results(solution_path)

    correct = verify(result, solution)
    return correct


if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_iev_1.txt'

    lines = readTextFile(path)
    genotype_pairs_count = [int(v) for v in lines[0].split()]
    expected = solve(genotype_pairs_count)

    print(expected)

    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, str(expected), 'w')

