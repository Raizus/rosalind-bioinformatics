

from typing import OrderedDict
import numpy as np
import numpy.typing as npt
from scipy.integrate import odeint

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.utils.utls import write_to_csv


def reaction_system(
    y,
    t: float,
    reactions: OrderedDict[int, Reaction],
    rate_constants: OrderedDict[int, float]
):
    # Create an array to store the rate of change for each species (initialized to zero)
    dydt = np.zeros_like(y)

    # Loop over each reaction
    for reaction_id, reaction in reactions.items():
        rate_constant = rate_constants[reaction_id]

        # Calculate the rate of the reaction based on the concentrations of reactants
        rate = rate_constant
        for reactant in reaction.reactants:
            # Multiply by the concentration of each reactant
            rate *= y[reactant]

        # Update the rate of change for reactants (decrease in concentration)
        for reactant in reaction.reactants:
            dydt[reactant] -= rate

        # Update the rate of change for products (increase in concentration)
        for product in reaction.products:
            dydt[product] += rate

    return dydt


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


class ODESimulator:
    cdat_filename: str = 'output.cdat'
    gdat_filename: str = 'output.gdat'

    def __init__(self,
                 cdat_filename: str = 'output.cdat',
                 gdat_filename: str = 'output.gdat'
                 ) -> None:
        self.cdat_filename = cdat_filename
        self.gdat_filename = gdat_filename

    def solve(self,
              concentrations: npt.NDArray[np.float_],
              t_span: npt.NDArray[np.float_],
              reactions: OrderedDict[int, Reaction],
              rate_constants: OrderedDict[int, float],
              groups: OrderedDict[int, ObservablesGroup]):

        sol = odeint(reaction_system, concentrations, t_span,
                     args=(reactions, rate_constants))

        n = concentrations.size
        # write cdat file
        header = ['Time'] + [f"S{i}" for i in range(n)]
        write_to_csv(self.cdat_filename, 'w', [header])
        aux = np.hstack((t_span.reshape((t_span.size, 1)), sol))
        write_to_csv(self.cdat_filename, 'a', aux)

        weights = get_groups_weight_matrix(groups, n)
        groups_conc = sol @ weights

        # write gdat file
        header = ['Time'] + [group.name for group in groups.values()]
        write_to_csv(self.gdat_filename, 'w', [header])
        aux = np.hstack((t_span.reshape((t_span.size, 1)), groups_conc))
        write_to_csv(self.gdat_filename, 'a', aux)
