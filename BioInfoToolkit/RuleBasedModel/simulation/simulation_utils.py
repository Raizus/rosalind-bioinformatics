from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup


import numpy as np


from typing import OrderedDict


def get_groups_weight_matrix(groups: OrderedDict[int, 'ObservablesGroup'], n: int):
    """
    Returns a n by m matrix (numpy array) where m is the length of groups, and 
    each column of the matrix is given by group.get_weights_array(n).

    Parameters:
    - groups: OrderedDict[int, ObservablesGroup] - a dictionary of groups
    - n: int - the length of the weights array for each group

    Returns:
    - A numpy array of shape (n, m)
    """
    # Number of groups (m) is the length of the dictionary
    m = len(groups)

    # Initialize an empty n x m matrix
    weight_matrix = np.empty((n, m), dtype=np.float64)

    # Iterate over groups and fill the matrix
    for idx, group in enumerate(groups.values()):
        weight_matrix[:, idx] = group.get_weights_array(n)

    return weight_matrix