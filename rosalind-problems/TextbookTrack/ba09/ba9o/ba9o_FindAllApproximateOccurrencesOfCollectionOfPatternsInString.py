from BioInfoToolkit.IO.IO import readTextFile, result_path_from_input_path, solution_path_from_input_path, writeTextFile
from BioInfoToolkit.Sequences.BarrowsWheeler import approximateMultiplePatternMatchingBwt
import os

from BioInfoToolkit.Sequences.StringUtils import hamming_distance


OutputT = list[int]


def verify(result: OutputT, solution: OutputT) -> bool:
    correct = result == solution
    return correct


def load_results(path: str) -> OutputT:
    lines = readTextFile(path)
    matches_idxs = [int(v) for v in lines[0].split()]
    return matches_idxs


def compare_result_solution(path_res: str, path_sol: str):
    result = load_results(path_res)
    solution = load_results(path_sol)

    correct = verify(result, solution)
    return correct



# FindAllApproximateOccurrencesOfCollectionOfPatternsInString
if __name__ == "__main__":
    cwd = os.path.realpath(os.path.dirname(__file__))
    path = f'{cwd}/rosalind_ba9o_1.txt'

    lines = readTextFile(path)
    string = lines[0] + '$'
    patterns = lines[1].split()
    d = int(lines[2])

    # matches = approximatePatternMatching(string, patterns[0], d)

    matches_dict = approximateMultiplePatternMatchingBwt(string, patterns, d)

    matches: list[int] = []
    for match_i in matches_dict.values():
        matches.extend(match_i)
    matches.sort()

    out = ' '.join(str(i) for i in matches)
    print(out)
    
    result_path = result_path_from_input_path(path)
    writeTextFile(result_path, out, 'w')

    # for i, pattern in enumerate(patterns):
    #     match_i = matches_dict[i]
    #     out_i = ' '.join(str(idx) for idx in match_i)
    #     print(f"{pattern}: {out_i}")
    #     writeTextFile(result_path, f"{pattern}: {out_i}", 'a')

    solution_path = solution_path_from_input_path(path)
    solution = load_results(solution_path)

    correct = verify(matches, solution)
    print(correct)

    # string = "panamabananas$"
    # pattern = "ana"
    # d = 1
    # approximatePatternMatching(string, pattern, d)
