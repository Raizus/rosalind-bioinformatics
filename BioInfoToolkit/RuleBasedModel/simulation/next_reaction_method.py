from typing import OrderedDict

import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction_generation import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import compute_propensities, \
    create_cdat, create_gdat, get_affected_reactions, get_groups_weight_matrix, get_species_to_reaction_map, update_affected_propensities, write_data_row
from BioInfoToolkit.RuleBasedModel.simulation.simulator import SimulatorABC


def update_taus(
    propensities: npt.NDArray[np.float_],
    new_propensities: npt.NDArray[np.float_],
    taus_es: npt.NDArray[np.float_],
    tau_min: float,
    reaction_idx: int
):
    for r_id, (old_rate, new_rate) in enumerate(zip(propensities, new_propensities)):
        new_rate = new_propensities[r_id]
        old_rate = propensities[r_id]
        if new_rate == 0:
            taus_es[r_id] = np.inf
        elif old_rate == 0 or r_id == reaction_idx:
            taus_es[r_id] = -np.log(np.random.rand()) / new_rate
        else:
            ratio = old_rate / new_rate
            taus_es[r_id] = ratio * (taus_es[r_id] - tau_min)

    return taus_es


class NextReactionMethod(SimulatorABC):

    def simulate(
        self,
        concentrations: OrderedDict[int, int],
        reactions: OrderedDict[int, Reaction],
        rate_constants: OrderedDict[int, float],
        groups: OrderedDict[int, ObservablesGroup]
    ):

        y = np.array(list(concentrations.values()), dtype=np.float64)
        t_start, y = self.get_initial(y)
        time = t_start

        num_species = len(concentrations)
        t_end = self.sim_params['t_end']
        n_steps = self.sim_params['n_steps']
        recording_times: npt.NDArray[np.float_] | None = None

        # If n_steps is provided, generate time points at which to record data
        if n_steps is not None:
            recording_times = np.linspace(t_start, t_end, n_steps)
            self.next_recording_idx = 1  # to track the next recording time index

        sp_to_reaction_map = get_species_to_reaction_map(reactions)
        weights = get_groups_weight_matrix(groups, num_species)

        # write cdat file
        create_cdat(self.cdat_filename, y.size)
        write_data_row(self.cdat_filename, time, y)

        # write to gdat file
        groups_conc = y @ weights
        create_gdat(self.gdat_filename, groups)
        write_data_row(self.gdat_filename, time, groups_conc)

        propensities = compute_propensities(y, reactions, rate_constants)
        taus_es = -np.log(np.random.rand(propensities.size)) / propensities

        while time < t_end:
            # Find the reaction with the smallest tau (next to happen)
            next_reaction_index = int(np.argmin(taus_es))
            tau_min = taus_es[next_reaction_index]

            # Advance time by the smallest tau
            time += tau_min

            # Apply the effect of the chosen reaction on concentrations
            reaction = reactions[next_reaction_index]
            species: set[int] = set()
            for reactant in reaction.reactants:
                y[reactant-1] -= 1
                species.add(reactant)
            for prod in reaction.products:
                y[prod-1] += 1
                species.add(prod)

            affected_rxns = get_affected_reactions(sp_to_reaction_map, species)

            # Recompute the rates after the reaction affects the system
            new_propensities = propensities.copy()
            new_propensities = update_affected_propensities(
                new_propensities, rate_constants,
                y, reactions, affected_rxns)

            # Update taus: only those reactions affected by the change need updates
            taus_es = update_taus(propensities, new_propensities, taus_es,
                                  tau_min, next_reaction_index)
            propensities = new_propensities

            self.record_line(time, y, weights, recording_times)

        return y
    