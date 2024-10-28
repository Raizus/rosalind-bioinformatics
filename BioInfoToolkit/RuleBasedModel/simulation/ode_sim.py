

from typing import OrderedDict
import numpy as np
import numpy.typing as npt
from scipy.integrate import odeint

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction_generation import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import create_cdat, create_gdat, \
    get_groups_weight_matrix
from BioInfoToolkit.RuleBasedModel.simulation.simulator import SimulatorABC
from BioInfoToolkit.RuleBasedModel.utils.utls import write_to_csv


def reaction_system(
    y: npt.NDArray[np.float_],
    t: float, # pylint: disable=unused-argument
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


class ODESimulator(SimulatorABC):

    def simulate(
        self,
        concentrations: npt.NDArray[np.float_],
        reactions: OrderedDict[int, Reaction],
        rate_constants: OrderedDict[int, float],
        groups: OrderedDict[int, ObservablesGroup]
    ):

        t_start = self.sim_params['t_start']
        t_end = self.sim_params['t_end']
        n_steps = self.sim_params['n_steps']
        n_steps = n_steps if n_steps is not None else 1000
        self.sim_params['n_steps'] = n_steps

        t_start, concentrations = self.get_initial(concentrations)
        t_span = np.linspace(t_start, t_end, n_steps)

        sol = odeint(
            reaction_system,
            concentrations,
            t_span,
            args=(reactions, rate_constants),
            rtol=self.sim_params['rtol'],
            atol=self.sim_params['atol']
        )

        n = concentrations.size

        # write cdat file
        create_cdat(self.cdat_filename, n)
        aux = np.hstack((t_span.reshape((t_span.size, 1)), sol))
        write_to_csv(self.cdat_filename, 'a', aux)

        weights = get_groups_weight_matrix(groups, n)
        groups_conc = sol @ weights

        # write gdat file
        create_gdat(self.gdat_filename, groups)
        aux = np.hstack((t_span.reshape((t_span.size, 1)), groups_conc))
        write_to_csv(self.gdat_filename, 'a', aux)

        y = sol[-1, :]
        return y
