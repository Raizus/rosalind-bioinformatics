

from typing import OrderedDict
import numpy as np
import numpy.typing as npt
from scipy.integrate import odeint

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import get_groups_weight_matrix
from BioInfoToolkit.RuleBasedModel.utils.utls import write_to_csv


def reaction_system(
    y: npt.NDArray[np.float_],
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
            rate *= y[reactant-1]

        # Update the rate of change for reactants (decrease in concentration)
        for reactant in reaction.reactants:
            dydt[reactant-1] -= rate

        # Update the rate of change for products (increase in concentration)
        for product in reaction.products:
            dydt[product-1] += rate

    return dydt


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
        header = ['Time'] + [f"S{i+1}" for i in range(n)]
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
