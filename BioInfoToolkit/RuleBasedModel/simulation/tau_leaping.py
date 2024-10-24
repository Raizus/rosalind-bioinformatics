from typing import OrderedDict

import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import calculate_tau_full, \
    compute_propensities, compute_stoichiometry_matrix, create_cdat, create_gdat, get_groups_weight_matrix, write_data_row
from BioInfoToolkit.RuleBasedModel.simulation.simulator import SimulatorABC


def update_concentrations(
    y: npt.NDArray[np.float_],
    reactions: OrderedDict[int, Reaction],
    reaction_counts: npt.NDArray[np.int_],
):
    # Update the concentrations based on how many times each reaction occurs
    for reaction_id, reaction in reactions.items():
        num_firings = reaction_counts[reaction_id]

        # Before applying the reaction, check if any reactants would go negative
        for reactant in reaction.reactants:
            if y[reactant-1] - num_firings < 0:
                # If the reaction would make the population negative, adjust num_firings
                num_firings = int(y[reactant-1])

        # Apply the adjusted number of firings to the concentrations
        for reactant in reaction.reactants:
            y[reactant-1] -= num_firings

        for product in reaction.products:
            y[product-1] += num_firings

    return y


class TauLeapingSimulator(SimulatorABC):

    def solve(self,
              concentrations: OrderedDict[int, int],
              reactions: OrderedDict[int, Reaction],
              rate_constants: OrderedDict[int, float],
              groups: OrderedDict[int, ObservablesGroup]):

        epsilon = 0.03

        num_species = len(concentrations)
        t_start = self.sim_params['t_start']
        t_end = self.sim_params['t_end']
        n_steps = self.sim_params['n_steps']
        given_tau = self.sim_params['tau']
        time = t_start
        recording_times: npt.NDArray[np.float_] | None = None

        y = np.array(list(concentrations.values()), dtype=np.float64)
        weights = get_groups_weight_matrix(groups, num_species)

        # If n_steps is provided, generate time points at which to record data
        if n_steps is not None:
            recording_times = np.linspace(t_start, t_end, n_steps)
            self.next_recording_idx = 1  # to track the next recording time index

        # write cdat file
        create_cdat(self.cdat_filename, y.size)
        write_data_row(self.cdat_filename, time, y)

        # write to gdat file
        groups_conc = y @ weights
        create_gdat(self.gdat_filename, groups)
        write_data_row(self.gdat_filename, time, groups_conc)

        stoichiometry_matrix = compute_stoichiometry_matrix(reactions, num_species)

        print("Starting simulation with tau-leaping algorithm...")
        print(f"\tt = {time}")

        # Perform tau-leaping simulation
        while time < t_end:
            # Calculate propensities for each reaction
            propensities = compute_propensities(y, reactions, rate_constants)

            # Calculate tau using the Cao method
            if given_tau is None:
                tau = calculate_tau_full(
                    y, reactions, propensities, stoichiometry_matrix, epsilon=epsilon)
            else:
                tau = given_tau

            # Draw the number of times each reaction occurs from a Poisson distribution
            reaction_counts = np.random.poisson(propensities * tau)

            # Update the concentrations based on how many times each reaction occurs
            y = update_concentrations(y, reactions, reaction_counts)

            # Update time
            time += tau

            self.record_line(time, y, weights, recording_times)

        return y
