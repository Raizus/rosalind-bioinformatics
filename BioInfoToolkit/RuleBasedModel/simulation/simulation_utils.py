from functools import reduce
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup


import numpy as np


from typing import OrderedDict

from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction


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


def compute_reaction_rates(
    rate_constants: OrderedDict[int, float],
    concentrations: OrderedDict[int, int],
    reactions: OrderedDict[int, Reaction]
):
    rates: OrderedDict[int, float] = OrderedDict()
    # compute the rates for each reaction
    # the rate is equal to the product of the concentrations of the left
    # side reagents and the reaction rate constant
    for r_id, rxn in reactions.items():
        rate_constant = rate_constants[r_id]
        concent_reactants: list[int] = []
        for reactant in rxn.reactants:
            if reactant in concentrations:
                concent_reactants.append(concentrations[reactant])
            else:
                raise KeyError(
                    f"Species with id '{reactant}' not in the concentrations dictionary.")
        rate = rate_constant * \
            reduce(lambda x, y: x*y, concent_reactants, 1.0)
        rates[r_id] = rate

    return rates