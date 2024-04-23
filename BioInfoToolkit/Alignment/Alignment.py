import csv
from itertools import combinations, product
import math
import sys
from typing import Callable
import drawsvg as draw
import numpy as np

PAM250: dict[tuple[str, str], int] = {
    ('A', 'A'): 2, ('A', 'C'): -2, ('A', 'D'): 0, ('A', 'E'): 0, ('A', 'F'): -3, ('A', 'G'): 1, ('A', 'H'): -1, ('A', 'I'): -1, ('A', 'K'): -1, ('A', 'L'): -2, ('A', 'M'): -1, ('A', 'N'): 0, ('A', 'P'): 1, ('A', 'Q'): 0, ('A', 'R'): -2, ('A', 'S'): 1, ('A', 'T'): 1, ('A', 'V'): 0, ('A', 'W'): -6, ('A', 'Y'): -3, ('C', 'A'): -2, ('C', 'C'): 12, ('C', 'D'): -5, ('C', 'E'): -5, ('C', 'F'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -2, ('C', 'K'): -5, ('C', 'L'): -6, ('C', 'M'): -5, ('C', 'N'): -4, ('C', 'P'): -3, ('C', 'Q'): -5, ('C', 'R'): -4, ('C', 'S'): 0, ('C', 'T'): -2, ('C', 'V'): -2, ('C', 'W'): -8, ('C', 'Y'): 0, ('D', 'A'): 0, ('D', 'C'): -5, ('D', 'D'): 4, ('D', 'E'): 3, ('D', 'F'): -6, ('D', 'G'): 1, ('D', 'H'): 1, ('D', 'I'): -2, ('D', 'K'): 0, ('D', 'L'): -4, ('D', 'M'): -3, ('D', 'N'): 2, ('D', 'P'): -1, ('D', 'Q'): 2, ('D', 'R'): -1, ('D', 'S'): 0, ('D', 'T'): 0, ('D', 'V'): -2, ('D', 'W'): -7, ('D', 'Y'): -4, ('E', 'A'): 0, ('E', 'C'): -5, ('E', 'D'): 3, ('E', 'E'): 4, ('E', 'F'): -5, ('E', 'G'): 0, ('E', 'H'): 1, ('E', 'I'): -2, ('E', 'K'): 0, ('E', 'L'): -3, ('E', 'M'): -2, ('E', 'N'): 1, ('E', 'P'): -1, ('E', 'Q'): 2, ('E', 'R'): -1, ('E', 'S'): 0, ('E', 'T'): 0, ('E', 'V'): -2, ('E', 'W'): -7, ('E', 'Y'): -4, ('F', 'A'): -3, ('F', 'C'): -4, ('F', 'D'): -6, ('F', 'E'): -5, ('F', 'F'): 9, ('F', 'G'): -5, ('F', 'H'): -2, ('F', 'I'): 1, ('F', 'K'): -5, ('F', 'L'): 2, ('F', 'M'): 0, ('F', 'N'): -3, ('F', 'P'): -5, ('F', 'Q'): -5, ('F', 'R'): -4, ('F', 'S'): -3, ('F', 'T'): -3, ('F', 'V'): -1, ('F', 'W'): 0, ('F', 'Y'): 7, ('G', 'A'): 1, ('G', 'C'): -3, ('G', 'D'): 1, ('G', 'E'): 0, ('G', 'F'): -5, ('G', 'G'): 5, ('G', 'H'): -2, ('G', 'I'): -3, ('G', 'K'): -2, ('G', 'L'): -4, ('G', 'M'): -3, ('G', 'N'): 0, ('G', 'P'): 0, ('G', 'Q'): -1, ('G', 'R'): -3, ('G', 'S'): 1, ('G', 'T'): 0, ('G', 'V'): -1, ('G', 'W'): -7, ('G', 'Y'): -5, ('H', 'A'): -1, ('H', 'C'): -3, ('H', 'D'): 1, ('H', 'E'): 1, ('H', 'F'): -2, ('H', 'G'): -2, ('H', 'H'): 6, ('H', 'I'): -2, ('H', 'K'): 0, ('H', 'L'): -2, ('H', 'M'): -2, ('H', 'N'): 2, ('H', 'P'): 0, ('H', 'Q'): 3, ('H', 'R'): 2, ('H', 'S'): -1, ('H', 'T'): -1, ('H', 'V'): -2, ('H', 'W'): -3, ('H', 'Y'): 0, ('I', 'A'): -1, ('I', 'C'): -2, ('I', 'D'): -2, ('I', 'E'): -2, ('I', 'F'): 1, ('I', 'G'): -3, ('I', 'H'): -2, ('I', 'I'): 5, ('I', 'K'): -2, ('I', 'L'): 2, ('I', 'M'): 2, ('I', 'N'): -2, ('I', 'P'): -2, ('I', 'Q'): -2, ('I', 'R'): -2, ('I', 'S'): -1, ('I', 'T'): 0, ('I', 'V'): 4, ('I', 'W'): -5, ('I', 'Y'): -1, ('K', 'A'): -1, ('K', 'C'): -5, ('K', 'D'): 0, ('K', 'E'): 0, ('K', 'F'): -5, ('K', 'G'): -2, ('K', 'H'): 0, ('K', 'I'): -2, ('K', 'K'): 5, ('K', 'L'): -3, ('K', 'M'): 0, ('K', 'N'): 1, ('K', 'P'): -1, ('K', 'Q'): 1, ('K', 'R'): 3, ('K', 'S'): 0, ('K', 'T'): 0, ('K', 'V'): -2, ('K', 'W'): -3, ('K', 'Y'): -4, ('L', 'A'): -2, ('L', 'C'): -6, ('L', 'D'): -4, ('L', 'E'): -3, ('L', 'F'): 2, ('L', 'G'): -4, ('L', 'H'): -2, ('L', 'I'): 2, ('L', 'K'): -3, ('L', 'L'): 6, ('L', 'M'): 4, ('L', 'N'): -3, ('L', 'P'): -3, ('L', 'Q'): -2, ('L', 'R'): -3, ('L', 'S'): -3, ('L', 'T'): -2, ('L', 'V'): 2, ('L', 'W'): -2, ('L', 'Y'): -1, ('M', 'A'): -1, ('M', 'C'): -5, ('M', 'D'): -3, ('M', 'E'): -2, ('M', 'F'): 0, ('M', 'G'): -3, ('M', 'H'): -2, ('M', 'I'): 2, ('M', 'K'): 0, ('M', 'L'): 4, ('M', 'M'): 6, ('M', 'N'): -2, ('M', 'P'): -2, ('M', 'Q'): -1, ('M', 'R'): 0, ('M', 'S'): -2, ('M', 'T'): -1, ('M', 'V'): 2, ('M', 'W'): -4, ('M', 'Y'): -2, ('N', 'A'): 0, ('N', 'C'): -4, ('N', 'D'): 2, ('N', 'E'): 1, ('N', 'F'): -3, ('N', 'G'): 0, ('N', 'H'): 2, ('N', 'I'): -2, ('N', 'K'): 1, ('N', 'L'): -3, ('N', 'M'): -2, ('N', 'N'): 2, ('N', 'P'): 0, ('N', 'Q'): 1, ('N', 'R'): 0, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'V'): -2, ('N', 'W'): -4, ('N', 'Y'): -2, ('P', 'A'): 1, ('P', 'C'): -3, ('P', 'D'): -1, ('P', 'E'): -1, ('P', 'F'): -5, ('P', 'G'): 0, ('P', 'H'): 0, ('P', 'I'): -2, ('P', 'K'): -1, ('P', 'L'): -3, ('P', 'M'): -2, ('P', 'N'): 0, ('P', 'P'): 6, ('P', 'Q'): 0, ('P', 'R'): 0, ('P', 'S'): 1, ('P', 'T'): 0, ('P', 'V'): -1, ('P', 'W'): -6, ('P', 'Y'): -5, ('Q', 'A'): 0, ('Q', 'C'): -5, ('Q', 'D'): 2, ('Q', 'E'): 2, ('Q', 'F'): -5, ('Q', 'G'): -1, ('Q', 'H'): 3, ('Q', 'I'): -2, ('Q', 'K'): 1, ('Q', 'L'): -2, ('Q', 'M'): -1, ('Q', 'N'): 1, ('Q', 'P'): 0, ('Q', 'Q'): 4, ('Q', 'R'): 1, ('Q', 'S'): -1, ('Q', 'T'): -1, ('Q', 'V'): -2, ('Q', 'W'): -5, ('Q', 'Y'): -4, ('R', 'A'): -2, ('R', 'C'): -4, ('R', 'D'): -1, ('R', 'E'): -1, ('R', 'F'): -4, ('R', 'G'): -3, ('R', 'H'): 2, ('R', 'I'): -2, ('R', 'K'): 3, ('R', 'L'): -3, ('R', 'M'): 0, ('R', 'N'): 0, ('R', 'P'): 0, ('R', 'Q'): 1, ('R', 'R'): 6, ('R', 'S'): 0, ('R', 'T'): -1, ('R', 'V'): -2, ('R', 'W'): 2, ('R', 'Y'): -4, ('S', 'A'): 1, ('S', 'C'): 0, ('S', 'D'): 0, ('S', 'E'): 0, ('S', 'F'): -3, ('S', 'G'): 1, ('S', 'H'): -1, ('S', 'I'): -1, ('S', 'K'): 0, ('S', 'L'): -3, ('S', 'M'): -2, ('S', 'N'): 1, ('S', 'P'): 1, ('S', 'Q'): -1, ('S', 'R'): 0, ('S', 'S'): 2, ('S', 'T'): 1, ('S', 'V'): -1, ('S', 'W'): -2, ('S', 'Y'): -3, ('T', 'A'): 1, ('T', 'C'): -2, ('T', 'D'): 0, ('T', 'E'): 0, ('T', 'F'): -3, ('T', 'G'): 0, ('T', 'H'): -1, ('T', 'I'): 0, ('T', 'K'): 0, ('T', 'L'): -2, ('T', 'M'): -1, ('T', 'N'): 0, ('T', 'P'): 0, ('T', 'Q'): -1, ('T', 'R'): -1, ('T', 'S'): 1, ('T', 'T'): 3, ('T', 'V'): 0, ('T', 'W'): -5, ('T', 'Y'): -3, ('V', 'A'): 0, ('V', 'C'): -2, ('V', 'D'): -2, ('V', 'E'): -2, ('V', 'F'): -1, ('V', 'G'): -1, ('V', 'H'): -2, ('V', 'I'): 4, ('V', 'K'): -2, ('V', 'L'): 2, ('V', 'M'): 2, ('V', 'N'): -2, ('V', 'P'): -1, ('V', 'Q'): -2, ('V', 'R'): -2, ('V', 'S'): -1, ('V', 'T'): 0, ('V', 'V'): 4, ('V', 'W'): -6, ('V', 'Y'): -2, ('W', 'A'): -6, ('W', 'C'): -8, ('W', 'D'): -7, ('W', 'E'): -7, ('W', 'F'): 0, ('W', 'G'): -7, ('W', 'H'): -3, ('W', 'I'): -5, ('W', 'K'): -3, ('W', 'L'): -2, ('W', 'M'): -4, ('W', 'N'): -4, ('W', 'P'): -6, ('W', 'Q'): -5, ('W', 'R'): 2, ('W', 'S'): -2, ('W', 'T'): -5, ('W', 'V'): -6, ('W', 'W'): 17, ('W', 'Y'): 0, ('Y', 'A'): -3, ('Y', 'C'): 0, ('Y', 'D'): -4, ('Y', 'E'): -4, ('Y', 'F'): 7, ('Y', 'G'): -5, ('Y', 'H'): 0, ('Y', 'I'): -1, ('Y', 'K'): -4, ('Y', 'L'): -1, ('Y', 'M'): -2, ('Y', 'N'): -2, ('Y', 'P'): -5, ('Y', 'Q'): -4, ('Y', 'R'): -4, ('Y', 'S'): -3, ('Y', 'T'): -3, ('Y', 'V'): -2, ('Y', 'W'): 0, ('Y', 'Y'): 10}

BLOSUM62: dict[tuple[str, str], int] = {
    ('A', 'A'): 4, ('A', 'C'): 0, ('A', 'D'): -2, ('A', 'E'): -1, ('A', 'F'): -2, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'K'): -1, ('A', 'L'): -1, ('A', 'M'): -1, ('A', 'N'): -2, ('A', 'P'): -1, ('A', 'Q'): -1, ('A', 'R'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'V'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('C', 'A'): 0, ('C', 'C'): 9, ('C', 'D'): -3, ('C', 'E'): -4, ('C', 'F'): -2, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'K'): -3, ('C', 'L'): -1, ('C', 'M'): -1, ('C', 'N'): -3, ('C', 'P'): -3, ('C', 'Q'): -3, ('C', 'R'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'V'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('D', 'A'): -2, ('D', 'C'): -3, ('D', 'D'): 6, ('D', 'E'): 2, ('D', 'F'): -3, ('D', 'G'): -1, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'K'): -1, ('D', 'L'): -4, ('D', 'M'): -3, ('D', 'N'): 1, ('D', 'P'): -1, ('D', 'Q'): 0, ('D', 'R'): -2, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'V'): -3, ('D', 'W'): -4, ('D', 'Y'): -3, ('E', 'A'): -1, ('E', 'C'): -4, ('E', 'D'): 2, ('E', 'E'): 5, ('E', 'F'): -3, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'K'): 1, ('E', 'L'): -3, ('E', 'M'): -2, ('E', 'N'): 0, ('E', 'P'): -1, ('E', 'Q'): 2, ('E', 'R'): 0, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'V'): -2, ('E', 'W'): -3, ('E', 'Y'): -2, ('F', 'A'): -2, ('F', 'C'): -2, ('F', 'D'): -3, ('F', 'E'): -3, ('F', 'F'): 6, ('F', 'G'): -3, ('F', 'H'): -1, ('F', 'I'): 0, ('F', 'K'): -3, ('F', 'L'): 0, ('F', 'M'): 0, ('F', 'N'): -3, ('F', 'P'): -4, ('F', 'Q'): -3, ('F', 'R'): -3, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'V'): -1, ('F', 'W'): 1, ('F', 'Y'): 3, ('G', 'A'): 0, ('G', 'C'): -3, ('G', 'D'): -1, ('G', 'E'): -2, ('G', 'F'): -3, ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'K'): -2, ('G', 'L'): -4, ('G', 'M'): -3, ('G', 'N'): 0, ('G', 'P'): -2, ('G', 'Q'): -2, ('G', 'R'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'V'): -3, ('G', 'W'): -2, ('G', 'Y'): -3, ('H', 'A'): -2, ('H', 'C'): -3, ('H', 'D'): -1, ('H', 'E'): 0, ('H', 'F'): -1, ('H', 'G'): -2, ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'K'): -1, ('H', 'L'): -3, ('H', 'M'): -2, ('H', 'N'): 1, ('H', 'P'): -2, ('H', 'Q'): 0, ('H', 'R'): 0, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'V'): -3, ('H', 'W'): -2, ('H', 'Y'): 2, ('I', 'A'): -1, ('I', 'C'): -1, ('I', 'D'): -3, ('I', 'E'): -3, ('I', 'F'): 0, ('I', 'G'): -4, ('I', 'H'): -3, ('I', 'I'): 4, ('I', 'K'): -3, ('I', 'L'): 2, ('I', 'M'): 1, ('I', 'N'): -3, ('I', 'P'): -3, ('I', 'Q'): -3, ('I', 'R'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'V'): 3, ('I', 'W'): -3, ('I', 'Y'): -1, ('K', 'A'): -1, ('K', 'C'): -3, ('K', 'D'): -1, ('K', 'E'): 1, ('K', 'F'): -3, ('K', 'G'): -2, ('K', 'H'): -1, ('K', 'I'): -3, ('K', 'K'): 5, ('K', 'L'): -2, ('K', 'M'): -1, ('K', 'N'): 0, ('K', 'P'): -1, ('K', 'Q'): 1, ('K', 'R'): 2, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'V'): -2, ('K', 'W'): -3, ('K', 'Y'): -2, ('L', 'A'): -1, ('L', 'C'): -1, ('L', 'D'): -4, ('L', 'E'): -3, ('L', 'F'): 0, ('L', 'G'): -4, ('L', 'H'): -3, ('L', 'I'): 2, ('L', 'K'): -2, ('L', 'L'): 4, ('L', 'M'): 2, ('L', 'N'): -3, ('L', 'P'): -3, ('L', 'Q'): -2, ('L', 'R'): -2, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'V'): 1, ('L', 'W'): -2, ('L', 'Y'): -1, ('M', 'A'): -1, ('M', 'C'): -1, ('M', 'D'): -3, ('M', 'E'): -2, ('M', 'F'): 0, ('M', 'G'): -3, ('M', 'H'): -2, ('M', 'I'): 1, ('M', 'K'): -1, ('M', 'L'): 2, ('M', 'M'): 5, ('M', 'N'): -2, ('M', 'P'): -2, ('M', 'Q'): 0, ('M', 'R'): -1, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'V'): 1, ('M', 'W'): -1, ('M', 'Y'): -1, ('N', 'A'): -2, ('N', 'C'): -3, ('N', 'D'): 1, ('N', 'E'): 0, ('N', 'F'): -3, ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'K'): 0, ('N', 'L'): -3, ('N', 'M'): -2, ('N', 'N'): 6, ('N', 'P'): -2, ('N', 'Q'): 0, ('N', 'R'): 0, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'V'): -3, ('N', 'W'): -4, ('N', 'Y'): -2, ('P', 'A'): -1, ('P', 'C'): -3, ('P', 'D'): -1, ('P', 'E'): -1, ('P', 'F'): -4, ('P', 'G'): -2, ('P', 'H'): -2, ('P', 'I'): -3, ('P', 'K'): -1, ('P', 'L'): -3, ('P', 'M'): -2, ('P', 'N'): -2, ('P', 'P'): 7, ('P', 'Q'): -1, ('P', 'R'): -2, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'V'): -2, ('P', 'W'): -4, ('P', 'Y'): -3, ('Q', 'A'): -1, ('Q', 'C'): -3, ('Q', 'D'): 0, ('Q', 'E'): 2, ('Q', 'F'): -3, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'K'): 1, ('Q', 'L'): -2, ('Q', 'M'): 0, ('Q', 'N'): 0, ('Q', 'P'): -1, ('Q', 'Q'): 5, ('Q', 'R'): 1, ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'V'): -2, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('R', 'A'): -1, ('R', 'C'): -3, ('R', 'D'): -2, ('R', 'E'): 0, ('R', 'F'): -3, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'K'): 2, ('R', 'L'): -2, ('R', 'M'): -1, ('R', 'N'): 0, ('R', 'P'): -2, ('R', 'Q'): 1, ('R', 'R'): 5, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'V'): -3, ('R', 'W'): -3, ('R', 'Y'): -2, ('S', 'A'): 1, ('S', 'C'): -1, ('S', 'D'): 0, ('S', 'E'): 0, ('S', 'F'): -2, ('S', 'G'): 0, ('S', 'H'): -1, ('S', 'I'): -2, ('S', 'K'): 0, ('S', 'L'): -2, ('S', 'M'): -1, ('S', 'N'): 1, ('S', 'P'): -1, ('S', 'Q'): 0, ('S', 'R'): -1, ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'V'): -2, ('S', 'W'): -3, ('S', 'Y'): -2, ('T', 'A'): 0, ('T', 'C'): -1, ('T', 'D'): -1, ('T', 'E'): -1, ('T', 'F'): -2, ('T', 'G'): -2, ('T', 'H'): -2, ('T', 'I'): -1, ('T', 'K'): -1, ('T', 'L'): -1, ('T', 'M'): -1, ('T', 'N'): 0, ('T', 'P'): -1, ('T', 'Q'): -1, ('T', 'R'): -1, ('T', 'S'): 1, ('T', 'T'): 5, ('T', 'V'): 0, ('T', 'W'): -2, ('T', 'Y'): -2, ('V', 'A'): 0, ('V', 'C'): -1, ('V', 'D'): -3, ('V', 'E'): -2, ('V', 'F'): -1, ('V', 'G'): -3, ('V', 'H'): -3, ('V', 'I'): 3, ('V', 'K'): -2, ('V', 'L'): 1, ('V', 'M'): 1, ('V', 'N'): -3, ('V', 'P'): -2, ('V', 'Q'): -2, ('V', 'R'): -3, ('V', 'S'): -2, ('V', 'T'): 0, ('V', 'V'): 4, ('V', 'W'): -3, ('V', 'Y'): -1, ('W', 'A'): -3, ('W', 'C'): -2, ('W', 'D'): -4, ('W', 'E'): -3, ('W', 'F'): 1, ('W', 'G'): -2, ('W', 'H'): -2, ('W', 'I'): -3, ('W', 'K'): -3, ('W', 'L'): -2, ('W', 'M'): -1, ('W', 'N'): -4, ('W', 'P'): -4, ('W', 'Q'): -2, ('W', 'R'): -3, ('W', 'S'): -3, ('W', 'T'): -2, ('W', 'V'): -3, ('W', 'W'): 11, ('W', 'Y'): 2, ('Y', 'A'): -2, ('Y', 'C'): -2, ('Y', 'D'): -3, ('Y', 'E'): -2, ('Y', 'F'): 3, ('Y', 'G'): -3, ('Y', 'H'): 2, ('Y', 'I'): -1, ('Y', 'K'): -2, ('Y', 'L'): -1, ('Y', 'M'): -1, ('Y', 'N'): -2, ('Y', 'P'): -3, ('Y', 'Q'): -1, ('Y', 'R'): -2, ('Y', 'S'): -2, ('Y', 'T'): -2, ('Y', 'V'): -1, ('Y', 'W'): 2, ('Y', 'Y'): 7}


def overlap_trim_s_t(s: str, t: str) -> tuple[str, str]:
    """Trims the prefix of s and t by the length of the inserts ('-') at the start of t,
    and trims the suffix of s and t by the length of the inserts ('-') at the end of s. Ex:
        s:      SSSSSXX-XX-X------
        t:      -----Y-YY-YYTTTTTT
        new_s:       XX-XX-X
        new_t:       Y-YY-YY
    Args:
        s (str):
        t (str):

    Returns:
        new_s, new_t (tuple[str, str]): 
    """
    m = len(s)
    new_s = s.rstrip('-')
    ds = m - len(new_s)
    t = t[:-ds]
    n = len(t)

    new_t = t.lstrip('-')
    dt = n - len(new_t)
    new_s = new_s[dt:]

    return new_s, new_t


class SimilarityScore:
    matchScore: int
    mismatchScore: int
    similarityDict: dict[tuple[str, str], int] | None

    def __init__(self, matchScore: int = 1, mismatchScore: int = -1,
                 similarityDict: dict[tuple[str, str], int] | None = None) -> None:
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore
        self.similarityDict = similarityDict

    def score(self, x: str, y: str) -> int:
        if self.similarityDict:
            score = self.similarityDict[(x, y)]
            return score

        score = self.matchScore if x == y else self.mismatchScore
        return score


defaultSimilarityScore = SimilarityScore()


def hamming_distance(seq1: str, seq2: str) -> int:
    """Computes hamming distance between two strings

    Args:
        seq1 (str): First string
        seq2 (str): Second string

    Returns:
        int: Hamming distance
    """
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def editDistance(seq1: str, seq2: str) -> int:
    """Computes the edit distance between two strings, using the Wagner-Fisher algorithm.

    https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm

    Args:
        seq1 (str):
        seq2 (str):

    Returns:
        int: dist
    """
    m = len(seq1)

    curr = [i for i in range(m+1)]

    for j, yj in enumerate(seq2):
        prev, curr = curr, [j+1]
        for i, xi in enumerate(seq1):
            curr.append(prev[i] if xi == yj else min(
                prev[i+1], prev[i], curr[-1]) + 1)

    return curr[-1]


def linearGapGlobalAlignmentCost(alignedSeq1: str, alignedSeq2: str, gapPenalty: int,
                                 similarityScore: SimilarityScore) -> int:
    """Computes the alignment score of two aligned sequences, using a linear gap penalty given by:
        gap(k) = k * gapPenalty, where k is the gap length (k>0)

    Args:
        alignedSeq1 (str): First Aligned Sequence
        alignedSeq2 (str): Second Aligned Sequence
        gapPenalty (int, optional): Linear gap penalty.
        similarityScore (SimilarityScore): 

    Returns:
        int: The alignment cost
    """

    cost = 0

    m = len(alignedSeq1)
    n = len(alignedSeq2)

    assert m == n, f"alignedSeq1 (m={m}) and alignedSeq2 (n={n}) must have the same length"

    for i, (xi, yi) in enumerate(zip(alignedSeq1, alignedSeq2)):
        if xi != '-' and yi != '-':
            match = similarityScore.score(xi, yi)
            cost += match
        elif xi == '-' and yi == '-':
            raise Exception(
                f"x[{i}] and y[{i}] are both gaps. Alignment is invalid")
        else:
            gapCost = gapPenalty
            cost += gapCost

    return cost


def constantGapGlobalAlignmentCost(alignedSeq1: str, alignedSeq2: str,
                                   gapPenalty: int, similarityScore: SimilarityScore) -> int:
    """Computes the alignment score of two aligned sequences, using a constant gap penalty. In a constant gap penalty, every gap receives some predetermined constant penalty, regardless of its length. Thus, the insertion or deletion of 1000 contiguous symbols is penalized equally to that of a single symbol.

    Args:
        alignedSeq1 (str): First Aligned Sequence
        alignedSeq2 (str): Second Aligned Sequence
        gapPenalty (int, optional): Constant gap penalty.
        similarityScore (SimilarityScore): 

    Returns:
        int: The alignment cost with constant gap penalty
    """
    cost = 0

    m = len(alignedSeq1)
    n = len(alignedSeq2)

    assert m == n, f"alignedSeq1 (m={m}) and alignedSeq2 (n={n}) must have the same length"

    x_gap_len = 0
    y_gap_len = 0
    for i, (xi, yi) in enumerate(zip(alignedSeq1, alignedSeq2)):
        if xi != '-' and yi != '-':
            match = similarityScore.score(xi, yi)
            cost += match
            x_gap_len = 0
            y_gap_len = 0
        elif xi == '-' and yi == '-':
            raise Exception(
                f"x[{i}] and y[{i}] are both gaps. Alignment is invalid")
        elif xi == '-':
            gapCost = gapPenalty if x_gap_len == 0 else 0
            cost += gapCost
            x_gap_len += 1
            y_gap_len = 0
        else:
            gapCost = gapPenalty if y_gap_len == 0 else 0
            cost += gapCost
            x_gap_len = 0
            y_gap_len += 1

    return cost


def affineGapGlobalAlignmentCost(alignedSeq1: str, alignedSeq2: str,
                                 similarityScore: SimilarityScore,
                                 gapOpeningPenalty: int,
                                 gapExtendingPenalty: int) -> int:
    """Computes the alignment score of two aligned sequences, using an affine gap penalty given by:
        gap(k) = gapOpeningPenalty + (k - 1) * gapExtendingPenalty, where k is the gap length (k>0)

    Args:
        alignedSeq1 (str): First Aligned Sequence
        alignedSeq2 (str): Second Aligned Sequence
        similarityScore (SimilarityScore): 
        gapOpeningPenalty (int, optional): Cost of opening a gap.
        gapExtendingPenalty (int, optional): Cost of extending a gap.

    Returns:
        int: The alignment cost with affine gap penalty
    """
    cost = 0
    m = len(alignedSeq1)
    n = len(alignedSeq2)

    assert m == n, f"alignedSeq1 (m={m}) and alignedSeq2 (n={n}) must have the same length"

    x_gap_len = 0
    y_gap_len = 0
    for i, (xi, yi) in enumerate(zip(alignedSeq1, alignedSeq2)):
        if xi != '-' and yi != '-':
            match = similarityScore.score(xi, yi)
            cost += match
            x_gap_len = 0
            y_gap_len = 0
        elif xi == '-':
            gapCost = gapOpeningPenalty if x_gap_len == 0 else gapExtendingPenalty
            cost += gapCost
            x_gap_len += 1
            y_gap_len = 0
        else:
            gapCost = gapOpeningPenalty if y_gap_len == 0 else gapExtendingPenalty
            cost += gapCost
            x_gap_len = 0
            y_gap_len += 1

    return cost


def linearGapSemiGlobalAlignmentCost(alignedSeq1: str, alignedSeq2: str, gapPenalty: int,
                                     similarityScore: SimilarityScore) -> int:
    """_summary_

    Args:
        alignedSeq1 (str): First Aligned Sequence
        alignedSeq2 (str): Second Aligned Sequence
        gapPenalty (int): _description_
        similarityScore (SimilarityScore): 

    Raises:
        Exception: _description_
        Exception: _description_
        Exception: _description_

    Returns:
        int: Alignment cost
    """
    cost = 0

    m = len(alignedSeq1)
    n = len(alignedSeq2)

    assert m == n, f"alignedSeq1 (m={m}) and alignedSeq2 (n={n}) must have the same length"

    # trim start
    ss1l = alignedSeq1.lstrip('-')
    ss2l = alignedSeq2.lstrip('-')
    d1 = m - len(ss1l)
    d2 = n - len(ss2l)
    if d1 != 0 and d2 != 0:
        raise Exception(
            "Sequences have aligned gaps at the start. Alignment is invalid")
    trimmed2 = ss2l[d1:]
    trimmed1 = ss1l[d2:]

    # trim end
    ss1r = trimmed1.rstrip('-')
    ss2r = trimmed2.rstrip('-')
    d1 = len(trimmed1) - len(ss1r)
    d2 = len(trimmed2) - len(ss2r)
    if d1 != 0 and d2 != 0:
        raise Exception(
            "Sequences have aligned gaps at the end. Alignment is invalid")
    trimmed2 = ss2r[:len(ss2r)-d1]
    trimmed1 = ss1r[:len(ss1r)-d2]

    for i, (xi, yi) in enumerate(zip(trimmed1, trimmed2)):
        if xi == '-' and yi == '-':
            raise Exception(
                f"x[{i}] and y[{i}] are both gaps. Alignment is invalid")
        elif xi != '-' and yi != '-':
            match = similarityScore.score(xi, yi)
            cost += match
        else:
            gapCost = gapPenalty
            cost += gapCost

    return cost


def linearGapOverlapAlignmentCost(s: str, t: str, gapPenalty: int, similarityScore: SimilarityScore) -> int | float:
    if s[-1] == '-' or t[0] == '-':
        return math.inf
    
    score = 0
    # new_s, new_t = overlap_trim_s_t(s,t)

    for i, (si, ti) in enumerate(zip(s, t)):
        if si == '-' and ti == '-':
            raise Exception(
                f"s[{i}] and t[{i}] are both gaps. Alignment is invalid")
        if si == '-' or ti == '-':
            score += gapPenalty
        else:
            match = similarityScore.score(si, ti)
            score += match
    return score


def globalAlignmentLinearGapGen(seq1: str, seq2: str, gap_penalty: int, similarity_score: SimilarityScore):
    m = len(seq1)
    n = len(seq2)

    last_line = [gap_penalty*j for j in range(n+1)]
    yield last_line

    for i, xi in enumerate(seq1, 1):
        curr_line = [0 for _ in range(m+1)]
        curr_line[0] = i*gap_penalty
        for j, yj in enumerate(seq2, 1):
        # scores of col
            score = similarity_score.score(xi, yj)

            match = last_line[j-1] + score
            delete = last_line[j] + gap_penalty
            insert = curr_line[j-1] + gap_penalty

            max_val = max(match, delete, insert)
            curr_line[j] = max_val

        yield curr_line
        last_line = curr_line


def globalAlignmentLinearGapPenalty(seq1: str, seq2: str, gapPenalty: int = -1,
                                    options: SimilarityScore = defaultSimilarityScore):
    """Needleman-Wunsch algorithm to compute optimal sequence alignment

    Args:
        seq1 (str): a string sequence
        seq2 (str): another string sequence

    Returns:
        tuple[list[list[int]], list[list[list[str]]]]: tuple with scoring matrix and backtracking matrix
    """
    # Match: The two letters at the current index are the same.
    # Mismatch: The two letters at the current index are different.
    # Indel (Insertion or Deletion): The best alignment involves one letter aligning to a gap in the other string.

    # default scoring system:
    # match: +1
    # mismatch/indel: -1

    m = len(seq1)
    n = len(seq2)

    # score matrix
    F = [[0 for _ in range(n+1)] for _ in range(m+1)]
    originMat: list[list[list[str]]] = [[[]
                                         for _ in range(n+1)] for _ in range(m+1)]

    for i in range(m+1):
        F[i][0] = gapPenalty*i
    for j in range(n+1):
        F[0][j] = gapPenalty*j

    for i, xi in enumerate(seq1, 1):
        for j, yj in enumerate(seq2, 1):
            similarityScore = options.score(xi, yj)

            match = F[i-1][j-1] + similarityScore
            delete = F[i-1][j] + gapPenalty
            insert = F[i][j-1] + gapPenalty

            aux = ['diag', 'up', 'left']
            vals = [match, delete, insert]
            max_val = max(vals)

            origin = [aux[l] for l, val in enumerate(vals) if val == max_val]
            originMat[i][j] = origin
            F[i][j] = max_val

    return F, originMat


def globalAlignmentLinearGapPenaltyBacktrack(seq1: str, seq2: str, F: list[list[int]],
                                             gapPenalty: int, similarityScore: SimilarityScore):
    """_summary_

    Args:
        seq1 (str): _description_
        seq2 (str): _description_
        F (list[list[int]]): _description_
        gapPenalty (int): _description_
        similarityScore (SimilarityScore): _description_

    Returns:
        _type_: _description_
    """
    alignedSeq1 = ""
    alignedSeq2 = ""

    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:
        xi = seq1[i-1]
        yj = seq2[j-1]
        if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + similarityScore.score(xi, yj):
            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            i, j = i - 1, j - 1

        elif i > 0 and F[i][j] == F[i-1][j] + gapPenalty:
            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1

        else:
            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            j = j - 1

    return alignedSeq1, alignedSeq2


def countOptimalGlobalAlignmentsLinearGap(seq1: str, seq2: str,
                                          gapPenalty: int,
                                          similarityScore: SimilarityScore,
                                          mod: int | None = None) -> int:
    sys.setrecursionlimit(1500)

    F, _ = globalAlignmentLinearGapPenalty(
        seq1, seq2, gapPenalty, similarityScore)
    cache: dict[tuple[int, int], int] = dict()

    def _countOptimalAlignments(i: int, j: int) -> int:

        if i == 0 or j == 0:
            return 1

        if cache.get((i, j)) != None:
            return cache[(i, j)]

        xi = seq1[i-1]
        yj = seq2[j-1]

        _count = 0
        if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + similarityScore.score(xi, yj):
            _count = _count + _countOptimalAlignments(i-1, j-1)
        if i > 0 and F[i][j] == F[i-1][j] + gapPenalty:
            _count = _count + _countOptimalAlignments(i-1, j)
        if j > 0 and F[i][j] == F[i][j-1] + gapPenalty:
            _count = _count + _countOptimalAlignments(i, j-1)

        _count = _count if mod is None else _count % mod
        cache[(i, j)] = _count

        return _count

    count = _countOptimalAlignments(len(seq1), len(seq2))
    return count


def globalAlignmentLinearGapPenaltyMiddleNode(seq1: str, seq2: str,
                                              gapPenalty: int,
                                              similarityScore: SimilarityScore):
    """Needleman-Wunsch algorithm to compute optimal sequence alignment

    Args:
        seq1 (str): a string sequence
        seq2 (str): another string sequence

    Returns:

    """
    # Match: The two letters at the current index are the same.
    # Mismatch: The two letters at the current index are different.
    # Indel (Insertion or Deletion): The best alignment involves one letter aligning to a gap in the other string.

    # default scoring system:
    # match: +1
    # mismatch/indel: -1

    # m x n matrix
    m = len(seq1)
    n = len(seq2)
    middle_col = math.floor(n/2)

    # return final_score, curr


def globalAlignmentLinearPenaltyGapLinearSpace(seq1: str, seq2: str,
                                               gapPenalty: int,
                                               similarityScore: SimilarityScore):

    def recurse(i1: int, i2: int, j1: int, j2: int):
        if j1 == j2:
            return
        if i1 == i2:
            return

        middleRow = math.floor((i2-i1)/2)
        costs_left = globalAlignmentLinearGapPenaltyScoreInLinearSpace(
            seq1[i1:middleRow], seq2[j1:j2], gapPenalty, similarityScore)
        costs_right =\
            globalAlignmentLinearGapPenaltyScoreInLinearSpace(
                seq1[middleRow-1:i2][::-1], seq2[j1:j2][::-1], gapPenalty, similarityScore)

        a = 0

    recurse(0, len(seq1), 0, len(seq2))


def globalAlignmentAffineGapPenalty(seq1: str, seq2: str, openingPenalty: int = -11,
                                    extendingPenalty: int = -1,
                                    options: SimilarityScore = defaultSimilarityScore):
    """Needleman-Wunsch algorithm to compute optimal global sequence alignment

    Args:
        seq1 (str): a string sequence
        seq2 (str): another string sequence

    Returns:
        tuple[str, str, list[list[int]]]: tuple with the augmented strings and scoring matrix
    """

    m = len(seq1)
    n = len(seq2)

    negInf = -10**10

    # score matrix
    # G[i][j] score of the best alignemnt of x[:i-1] and y[:j-1] ending with character match or mismatch
    # X[i][j] score of the best alignemnt of x[:i-1] and y[:j-1] ending with space in X
    # Y[i][j] score of the best alignemnt of x[:i-1] and y[:j-1] ending with space in Y
    V = [[0 for _ in range(n+1)] for _ in range(m+1)]
    X = [[0 for _ in range(n+1)] for _ in range(m+1)]
    Y = [[0 for _ in range(n+1)] for _ in range(m+1)]
    G = [[0 for _ in range(n+1)] for _ in range(m+1)]

    for i in range(1, m+1):
        V[i][0] = openingPenalty + (i-1) * extendingPenalty
        X[i][0] = V[i][0]
        Y[i][0] = negInf
    for j in range(1, n+1):
        V[0][j] = openingPenalty + (j-1) * extendingPenalty
        Y[0][j] = V[0][j]
        X[0][j] = negInf
    for i, xi in enumerate(seq1, 1):
        for j, yj in enumerate(seq2, 1):
            match = options.score(xi, yj)

            X[i][j] = max(X[i][j-1] + extendingPenalty, V[i][j-1] +
                          openingPenalty)
            Y[i][j] = max(Y[i-1][j] + extendingPenalty,
                          V[i-1][j] + openingPenalty)
            G[i][j] = match + V[i-1][j-1]

            V[i][j] = max(G[i][j], X[i][j], Y[i][j])

    score = V[m][n]

    alignedSeq1 = ""
    alignedSeq2 = ""

    i = len(seq1)
    j = len(seq2)

    prev = 'diagonal'
    if V[i][j] == Y[i][j]:
        prev = 'up'
    if V[i][j] == X[i][j]:
        prev = 'left'

    while i > 0 or j > 0:
        xi = seq1[i-1]
        yj = seq2[j-1]

        up = max(Y[i-1][j] + extendingPenalty,
                 V[i-1][j] + openingPenalty)
        left = max(X[i][j-1] + extendingPenalty,
                   V[i][j-1] + openingPenalty)

        if prev == 'up' or j == 0:
            if Y[i][j] == V[i-1][j] + openingPenalty:
                prev = 'diagonal'
            else: 
                prev = 'up'

            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1
        elif prev == 'left' or i == 0:
            if X[i][j] == V[i][j-1] + openingPenalty:
                prev = 'diagonal'
            else:
                prev = 'left'

            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            j = j - 1
        else:
            if up == V[i][j]:
                prev = 'up'
            elif left > G[i][j]:
                prev = 'left'
            else:
                alignedSeq1 = xi + alignedSeq1
                alignedSeq2 = yj + alignedSeq2
                i, j = i - 1, j - 1

    return alignedSeq1, alignedSeq2, score, V, G, X, Y


def globalAlignmentConstantGapPenalty(seq1: str, seq2: str, gapPenalty: int,
                                      similarityScore: SimilarityScore):
    """Needleman-Wunsch algorithm to compute optimal sequence alignment, using a constant gap penalty.
    In a constant gap penalty, every gap receives some predetermined constant penalty, regardless of its length.

    http://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
    https://biochem218.stanford.edu/Projects%202004/Chan.pdf

    Args:
        seq1 (str): _description_
        seq2 (str): _description_
        gapPenalty (int, optional): constant gap penalty.
        similarityScore (SimilarityScore, optional):

    Returns:
        _type_: _description_
    """

    m = len(seq1)
    n = len(seq2)

    # score matrix
    # M[i][j] score of the best alignemnt of x[:i-1] and y[:j-1] ending with character match or mismatch
    M = [[0 for _ in range(n+1)] for _ in range(m+1)]
    # keeps track of the maximum in the column j
    max_j = [M[0][j] for j in range(1, n+1)]

    for i in range(1, m+1):
        M[i][0] = gapPenalty
    for j in range(1, n+1):
        M[0][j] = gapPenalty
    for i, xi in enumerate(seq1, 1):
        max_i = M[i][0]  # maximum in the current row
        for j, yj in enumerate(seq2, 1):
            match = similarityScore.score(xi, yj)

            a = match + M[i-1][j-1]
            b1 = max(max_i, M[i][j-1])
            b = b1 + gapPenalty
            # b = max(M[i][j2] for j2 in range(j)) + gapPenalty
            max_i = b1

            c1 = max(max_j[j-1], M[i-1][j])
            c = c1 + gapPenalty
            max_j[j-1] = c1
            # c = max(M[i2][j] for i2 in range(i)) + gapPenalty

            M[i][j] = max(a, b, c)

    score = M[m][n]
    return score, M


def alignmentMatrixToCsv(seq1: str, seq2: str, mat: list[list[int]], csv_path: str):
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(list('  ' + seq2))
        seq1_aux = ' ' + seq1
        for i, row in enumerate(mat):
            writer.writerow([seq1_aux[i]]+row)


def localAlignmentLinearGap(seq1: str, seq2: str, gapPenalty: int,
                            similarityScore: SimilarityScore = defaultSimilarityScore):
    """Performs a local alignment of two strings using the smith-waterman algorithm

    Args:
        seq1 (str): _description_
        seq2 (str): _description_
        gapPenalty (int): _description_
        similarityScore (SimilarityScore, optional): _description_. Defaults to defaultSimilarityScore.

    Returns:
        _type_: _description_
    """
    m = len(seq1)
    n = len(seq2)

    # score matrix
    H = [[0 for _ in range(n+1)] for _ in range(m+1)]
    max_val = 0
    max_ind: list[tuple[int, int]] = []

    for i, xi in enumerate(seq1, 1):
        for j, yj in enumerate(seq2, 1):
            score = similarityScore.score(xi, yj)

            match = H[i-1][j-1] + score
            v1 = H[i-1][j] + gapPenalty
            v2 = H[i][j-1] + gapPenalty
            # v1 = max(H[i-k][j] + gapPenalty * k for k in range(1,i+1))
            # v2 = max(H[i][j-l] + gapPenalty * l for l in range(1,j+1))

            H[i][j] = max(match, v1, v2, 0)

            if H[i][j] > max_val:
                max_val = H[i][j]
                max_ind = [(i, j)]
            elif H[i][j] == max_val:
                max_ind.append((i, j))

    # traceback
    i, j = max_ind[0]

    alignedSeq1 = ''
    alignedSeq2 = ''

    while H[i][j] > 0:
        xi = seq1[i-1]
        yj = seq2[j-1]

        if H[i][j] == H[i-1][j] + gapPenalty:
            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1
        elif H[i][j] == H[i][j-1] + gapPenalty:
            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            j = j - 1
        else:
            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            i, j = i - 1, j - 1

    return alignedSeq1, alignedSeq2, H, max_ind


def smithWatermanAffineGap(seq1: str, seq2: str, openingPenalty: int = -11,
                           extendingPenalty: int = -1, similarityScore: SimilarityScore = defaultSimilarityScore):
    m = len(seq1)
    n = len(seq2)

    # score matrix
    # G[i][j] score of the best alignemnt of x[:i-1] and y[:j-1] ending with character match or mismatch
    # X[i][j] score of the best alignemnt of x[:i-1] and y[:j-1] ending with space in X
    # Y[i][j] score of the best alignemnt of x[:i-1] and y[:j-1] ending with space in Y
    V = [[0 for _ in range(n+1)] for _ in range(m+1)]
    X = [[0 for _ in range(n+1)] for _ in range(m+1)]
    Y = [[0 for _ in range(n+1)] for _ in range(m+1)]

    max_val = 0
    max_ind: list[tuple[int, int]] = []

    for i, xi in enumerate(seq1, 1):
        for j, yj in enumerate(seq2, 1):
            match = similarityScore.score(xi, yj)

            X[i][j] = max(X[i][j-1] + extendingPenalty, V[i][j-1] +
                          openingPenalty)
            Y[i][j] = max(Y[i-1][j] + extendingPenalty,
                          V[i-1][j] + openingPenalty)
            Gij = match + V[i-1][j-1]

            V[i][j] = max(Gij, X[i][j], Y[i][j], 0)

            if V[i][j] > max_val:
                max_val = V[i][j]
                max_ind = [(i, j)]
            elif V[i][j] == max_val:
                max_ind.append((i, j))

    print(max_val)
    print(max_ind)

    # traceback
    i, j = max_ind[0]

    alignedSeq1 = ''
    alignedSeq2 = ''

    prev = 'diagonal'
    if V[i][j] == Y[i][j]:
        prev = 'up'
    if V[i][j] == X[i][j]:
        prev = 'left'

    while V[i][j] > 0:
        xi = seq1[i-1]
        yj = seq2[j-1]

        up = max(Y[i-1][j] + extendingPenalty,
                 V[i-1][j] + openingPenalty)
        left = max(X[i][j-1] + extendingPenalty,
                   V[i][j-1] + openingPenalty)

        if prev == 'up' or j == 0:
            if Y[i][j] == V[i-1][j] + openingPenalty:
                prev = 'diagonal'
            # else still 'up'

            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1
        elif prev == 'left' or i == 0:
            if X[i][j] == V[i][j-1] + openingPenalty:
                prev = 'diagonal'
            # else still 'left'

            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            j = j - 1
        else:
            match = similarityScore.score(xi, yj)
            if up == V[i][j]:
                prev = 'up'
            elif left > match + V[i-1][j-1]:
                prev = 'left'
            else:
                alignedSeq1 = xi + alignedSeq1
                alignedSeq2 = yj + alignedSeq2
                i, j = i - 1, j - 1

    del X, Y
    return alignedSeq1, alignedSeq2, V, max_ind


def fittingAlignmentLinearGap(seq: str, motif: str, gapPenalty: int, 
                              similarityScore: SimilarityScore):
    """Computes the fitting aligment between a string s and a motif t. A fitting alignment is an alignment of a substring of s against all of t.

    Args:
        seq (str): _description_
        motif (str): _description_
        gapPenalty (int): _description_
        similarityScore (SimilarityScore, optional): _description_. Defaults to defaultSimilarityScore.
    """
    m = len(seq)
    n = len(motif)

    H = [[0 for _ in range(n+1)] for _ in range(m+1)]
    originMat: list[list[list[str]]] = [[[] for _ in range(n+1)] for _ in range(m+1)]

    for j in range(1, n+1):
        H[0][j] = gapPenalty * j
        originMat[0][j] = ['left']

    aux = ['diag', 'up', 'left']
    for i, si in enumerate(seq, 1):
        for j, tj in enumerate(motif, 1):
            score = similarityScore.score(si, tj)

            match = H[i-1][j-1] + score
            delete = H[i-1][j] + gapPenalty
            insert = H[i][j-1] + gapPenalty

            vals = [match, delete, insert]
            H[i][j] = max(match, delete, insert)
            
            max_val = max(vals)
            origin = [aux[l] for l, val in enumerate(vals) if val == max_val]
            originMat[i][j] = origin


    max_score = max(H[i][n] for i in range(m+1))
    max_idxs = [i for i in range(m+1) if H[i][n] == max_score]

    return H, originMat, max_score


def fittingAlignmentLinearGapBacktrack(seq: str, motif: str, H: list[list[int]],
                                       gapPenalty: int, 
                                       similarityScore: SimilarityScore):
    """Computes the fitting aligment between a string s and a motif t. A fitting alignment is an alignment of a substring of s against all of t.

    Args:
        seq (str): _description_
        motif (str): _description_
        H (list[list[int]]): alignment cost matrix
        gapPenalty (int): _description_
        similarityScore (SimilarityScore, optional): _description_. Defaults to defaultSimilarityScore.
    """
    m = len(seq)
    n = len(motif)

    max_score = max(H[i][n] for i in range(m+1))
    max_idxs = [i for i in range(m+1) if H[i][n] == max_score]

    # backtrack
    alignedSeq1: str = ''
    alignedSeq2: str = ''

    i = max_idxs[0]
    j = n

    while j > 0:
        si = seq[i-1]
        tj = motif[j-1]

        if H[i][j] == H[i-1][j] + gapPenalty:
            alignedSeq1 = si + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1
        elif H[i][j] == H[i][j-1] + gapPenalty:
            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = tj + alignedSeq2
            j = j - 1
        else:
            alignedSeq1 = si + alignedSeq1
            alignedSeq2 = tj + alignedSeq2
            i, j = i - 1, j - 1

    return alignedSeq1, alignedSeq2


def semiglobalAligmentLinearGap(seq1: str, seq2: str,
                                gapPenalty: int,
                                similarityScore: SimilarityScore = defaultSimilarityScore):
    """A semiglobal alignment of strings s and t is an alignment in which any gaps appearing as prefixes or suffixes of s and t do not contribute to the alignment score.

    Args:
        seq1 (str): _description_
        seq2 (str): _description_
        gapPenalty (int): _description_
        similarityScore (SimilarityScore, optional): _description_. Defaults to defaultSimilarityScore.

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    """
    m = len(seq1)
    n = len(seq2)

    H = [[0 for _ in range(n+1)] for _ in range(m+1)]

    for i, xi in enumerate(seq1, 1):
        for j, yj in enumerate(seq2, 1):
            score = similarityScore.score(xi, yj)

            match = H[i-1][j-1] + score
            v1 = H[i-1][j] + (0 if j == n else gapPenalty)
            v2 = H[i][j-1] + (0 if i == m else gapPenalty)

            H[i][j] = max(match, v1, v2)

    score = H[m][n]
    i, j = m, n

    # backtrack
    alignedSeq1 = ''
    alignedSeq2 = ''

    while i > 0 or j > 0:
        xi = seq1[i-1]
        yj = seq2[j-1]

        if i > 0 and H[i][j] == H[i-1][j] + (0 if j == n or j == 0 else gapPenalty):
            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1
        elif j > 0 and H[i][j] == H[i][j-1] + (0 if i == m or i == 0 else gapPenalty):
            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            j = j - 1
        elif H[i][j] == H[i-1][j-1] + similarityScore.score(xi, yj):
            alignedSeq1 = xi + alignedSeq1
            alignedSeq2 = yj + alignedSeq2
            i, j = i - 1, j - 1
        else:
            raise Exception("Something went wrong.")

    return alignedSeq1, alignedSeq2, H, score


def overlapAlignmentLinearGap(s: str, t: str,
                              gapPenalty: int,
                              similarityScore: SimilarityScore):
    """An overlap alignment between two strings s and t is a local alignment of a suffix of s with a prefix of t. An optimal overlap alignment will therefore maximize an alignment score over all such substrings of s and t.

    Args:
        seq1 (str): _description_
        seq2 (str): _description_
        gapPenalty (int): _description_
        similarityScore (SimilarityScore):
    """
    m = len(s)
    n = len(t)

    H = [[0 for _ in range(n+1)] for _ in range(m+1)]

    for i in range(1, m):
        # gaps at the start of s are free
        H[i][0] = 0
    for j in range(1, n):
        H[0][j] = gapPenalty * j

    for i, si in enumerate(s, 1):
        for j, tj in enumerate(t, 1):
            score = similarityScore.score(si, tj)

            match = H[i-1][j-1] + score
            # score of best alignment of s[1..i] and t[1..j] ending with a space in s
            v1 = H[i-1][j] + gapPenalty
            # score of best alignment of s[1..i] and t[1..j] ending with a space in t
            # want to align suffix of s with prefix of t so gaps at the end of t are free
            v2 = H[i][j-1] + (0 if i == m else gapPenalty)

            H[i][j] = max(match, v1, v2)

    # last_line = [(gapPenalty*j, "", "-"*j) for j in range(n+1)]
    # for i, si in enumerate(s, 1):
    #     curr_line = [(0, "", "")]
    #     for j, tj in enumerate(t, 1):
    #         score = similarityScore.score(si, tj)
            
    #         match = last_line[j-1][0] + score
    #         v1 = last_line[j][0] + gapPenalty
    #         v2 = curr_line[j-1][0] + (0 if i == m else gapPenalty)


    score = H[m][n]
    i, j = m, n

    # backtrack
    aligned_s = ''
    aligned_t = ''

    while i > 0 or j > 0:
        si = s[i-1]
        tj = t[j-1]

        if i > 0 and H[i][j] == H[i-1][j] + (0 if j == 0 else gapPenalty):
            aligned_s = si + aligned_s
            aligned_t = '-' + aligned_t
            i = i - 1
        elif j > 0 and H[i][j] == H[i][j-1] + (0 if i == m else gapPenalty):
            aligned_s = '-' + aligned_s
            aligned_t = tj + aligned_t
            j = j - 1
        elif H[i][j] == H[i-1][j-1] + similarityScore.score(si, tj):
            aligned_s = si + aligned_s
            aligned_t = tj + aligned_t
            i, j = i - 1, j - 1
        else:
            raise Exception("Something went wrong.")

    return aligned_s, aligned_t, H, score


def drawAlignmentMatrixSvg(F: list[list[int]], seq1: str, seq2: str):
    """Draws the scoring matrix of an alignment onto an svg

    Args:
        F (list[list[int]]): scoring matrix
        seq1 (str): sequence 1
        seq2 (str): sequence 2

    Returns:
        (Drawing): the svg object
    """
    m = len(F)
    n = len(F[0])

    cell_size = 50
    font_size = 25
    height = cell_size*(m+1)
    width = cell_size*(n+1)
    d = draw.Drawing(width, height)
    d.append(draw.Rectangle(0, 0, width, height, fill="white"))

    # draw table lines
    for i in range(m+2):
        d.append(draw.Line(0, i*cell_size, width, i *
                 cell_size, stroke_width=2, stroke="black"))
    for j in range(n+2):
        d.append(draw.Line(j*cell_size, 0, j*cell_size,
                 height, stroke_width=2, stroke="black"))

    # draw sequence labels
    for i, c in enumerate(' '+seq1):
        d.append(draw.Text(c, font_size, (0.5)*cell_size, (i+1+0.5)
                 * cell_size, text_anchor='middle', center=True))
    for j, c in enumerate(' '+seq2):
        d.append(draw.Text(c, font_size, (j+1+0.5)*cell_size, (0.5)
                 * cell_size, text_anchor='middle', center=True))

    # draw matrix values
    for i, line in enumerate(F):
        for j, val in enumerate(line):
            d.append(draw.Text(str(val), font_size, (j+1+0.5)*cell_size,
                     (i+1+0.5)*cell_size, text_anchor='middle', center=True))

    return d


def add_backtrack_paths_to_alignment_svg_global_alignement(d: draw.Drawing,
                                                           originMat: list[list[list[str]]],
                                                           seq1: str, seq2: str, idxs: list[tuple[int,int]] = []):

    cell_size = 50
    if len(idxs) == 0:
        i = len(seq1)
        j = len(seq2)
        idxs.append((i,j))

    def dir_to_cords(i: int, j: int, direction: str):
        if direction == 'up':
            return (i-1, j)
        if direction == 'left':
            return (i, j-1)
        return (i-1, j-1)

    drawn_dict: dict[tuple[int, int], bool] = dict()

    def draw_backtrack_arrows(i: int, j: int):
        while i > 0 or j > 0:
            if drawn_dict.get((i, j), False):
                return
            for direction in originMat[i][j]:
                i2, j2 = dir_to_cords(i, j, direction)

                d.append(draw.Line((j+1+0.5)*cell_size, (i+1+0.5)*cell_size, (j2+1+0.5)
                         * cell_size, (i2+1+0.5)*cell_size, stroke_width=2, stroke="black"))

                draw_backtrack_arrows(i2, j2)
            drawn_dict[(i, j)] = True

    for i0, j0 in idxs:
        draw_backtrack_arrows(i0, j0)


def globalAlignmentLinearGapPenaltyScoreInLinearSpace2(seq1: str, seq2: str,
                                                       gapPenalty: int,
                                                       similarityScore: SimilarityScore):
    """Needleman-Wunsch algorithm to compute optimal sequence alignment

    Args:
        seq1 (str): a string sequence
        seq2 (str): another string sequence

    Returns:

    """
    # Match: The two letters at the current index are the same.
    # Mismatch: The two letters at the current index are different.
    # Indel (Insertion or Deletion): The best alignment involves one letter aligning to a gap in the other string.

    # default scoring system:
    # match: +1
    # mismatch/indel: -1

    # m x n matrix
    m = len(seq1)
    n = len(seq2)

    # scores of col
    prev = [0 for _ in range(m+1)]

    for i in range(m+1):
        prev[i] = gapPenalty*i

    for j, yj in enumerate(seq2, 1):
        # scores of col
        curr = [0 for _ in range(m+1)]
        curr[0] = j*gapPenalty
        for i, xi in enumerate(seq1, 1):
            score = similarityScore.score(xi, yj)

            match = prev[i-1] + score
            delete = prev[i] + gapPenalty
            insert = curr[i-1] + gapPenalty

            max_val = max(match, delete, insert)

            curr[i] = max_val

        prev = curr

    # i, j = m, n
    # while j == n:
    #     if i > 0 and curr[i] == curr[i-1] + gapPenalty:
    #         i = i - 1
    #     elif i > 0 and j > 0 and curr[i] == prev[i-1] + similarityScore.score(seq1[i-1], seq2[j-1]):
    #         i, j = i - 1, j - 1
    #     else:
    #         j = j - 1
    # middle_edge = (i, j)

    return prev



def middle_edge(seq1: str, seq2: str, gapPenalty: int,
                similarityScore: SimilarityScore):

    mid_col = len(seq2) // 2

    def meet(from_source_cost: list[int], to_sink_cost: list[int]):
        pair = (0, 0)
        best_score = -math.inf
        try:
            for i, value in enumerate(from_source_cost):
                horizontal = value + to_sink_cost[i] + gapPenalty
                if horizontal >= best_score:
                    best_score = horizontal
                    pair = (i, i)
                diagonal = value + \
                    to_sink_cost[i+1] + \
                    similarityScore.score(seq1[i], seq2[mid_col])
                if diagonal >= best_score:
                    best_score = diagonal
                    pair = (i, i+1)

        except IndexError:
            pass

        return pair, int(best_score)

    prefix = seq2[:mid_col]
    suffix = seq2[mid_col+1:]
    from_source_cost = globalAlignmentLinearGapPenaltyScoreInLinearSpace2(
        seq1, prefix, gapPenalty, similarityScore)

    seq1_b = seq1[::-1]
    seq2_b = suffix[::-1]
    to_sink_cost = globalAlignmentLinearGapPenaltyScoreInLinearSpace2(
        seq1_b, seq2_b, gapPenalty, similarityScore)

    to_sink_cost = to_sink_cost[::-1]

    pair, best_score = meet(from_source_cost, to_sink_cost)
    edge = ((pair[0], mid_col), (pair[1], mid_col+1))

    return edge, best_score


def globalAlignmentLinearGapPenaltyScoreInLinearSpace(seq1: str, seq2: str,
                                                      gapPenalty: int,
                                                      similarityScore: SimilarityScore):
    edges = []

    def recurse(top: int, bottom: int, left: int, right: int):
        if left == right:
            return 0
        if top == bottom:
            return 0

        edge, score = middle_edge(seq1[top:bottom], seq2[left:right], gapPenalty, similarityScore)

        n1 = (top + edge[0][0], left + edge[0][1])
        n2 = (top + edge[1][0], left + edge[1][1])
        edges.append((n1, n2))

        recurse(top, n1[0], left, n1[1])
        recurse(n2[0], bottom, n2[1], right)

        return score

    score = recurse(0, len(seq1), 0, len(seq2))

    edges.sort(key=lambda x: x[0][1])
    edges.sort(key=lambda x: x[0][0])

    m = len(seq1)
    n = len(seq2)
    nodes: list[tuple[int, int]] = []
    for edge in edges:
        n1 = edge[0]
        n2 = edge[1]

        if len(nodes) == 0 or nodes[-1] != n1:
            nodes.append(n1)
        nodes.append(n2)

    if nodes[0] != (0, 0):
        nodes = [(0, 0)] + nodes

    if nodes[-1] != (m, n):
        nodes.append((m, n))

    def alignmentFromEdges():
        alignedSeq1 = ''
        alignedSeq2 = ''
        idx = len(nodes)-1
        i, j = nodes[idx]
        while i > 0 or j > 0:
            (i2, j2) = nodes[idx-1]

            if i-i2 == 0:  # horizontal
                alignedSeq1 = '-'*(j-j2) + alignedSeq1
                alignedSeq2 = seq2[j2:j] + alignedSeq2
            elif j-j2 == 0:  # vertical
                alignedSeq1 = seq1[i2:i] + alignedSeq1
                alignedSeq2 = '-'*(i-i2) + alignedSeq2
            else:  # diagonal
                alignedSeq1 = seq1[i2:i] + alignedSeq1
                alignedSeq2 = seq2[j2:j] + alignedSeq2

            i, j = i2, j2
            idx -= 1

        return alignedSeq1, alignedSeq2

    alignedSeq1, alignedSeq2 = alignmentFromEdges()
    return score, alignedSeq1, alignedSeq2
