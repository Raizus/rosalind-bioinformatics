
def float_matrices_match(m1: list[list[float]], m2: list[list[float]]) -> bool:
    if len(m1) != len(m2):
        return False
    for line1, line2 in zip(m1, m2):
        if len(line1) != len(line2) or not all(abs(v1-v2) <= 0.001 for v1, v2 in zip(line1, line2)):
            return False
    return True
