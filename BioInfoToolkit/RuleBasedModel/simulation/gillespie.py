
from typing import OrderedDict
import random
import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import compute_propensities, \
    create_cdat, create_gdat, get_affected_reactions, get_groups_weight_matrix, get_species_to_reaction_map, update_affected_propensities, write_data_row
from BioInfoToolkit.RuleBasedModel.simulation.simulator import SimulatorABC


class GillespieSimulator(SimulatorABC):

    def simulate(self,
                 reactions: OrderedDict[int, Reaction],
                 rate_constants: OrderedDict[int, float],
                 concentrations: OrderedDict[int, int],
                 groups: OrderedDict[int, ObservablesGroup]):

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

        # Initialize observable groups concentrations
        weights = get_groups_weight_matrix(groups, num_species)

        # Write the header and first row of .cdat file
        create_cdat(self.cdat_filename, num_species)
        write_data_row(self.cdat_filename, time, y)

        # write to gdat file
        groups_conc = y @ weights
        create_gdat(self.gdat_filename, groups)
        write_data_row(self.gdat_filename, time, groups_conc)

        print("Starting simulation with Gillespie algorithm...")
        print(f"\tt = {time}")

        r_keys = list(reactions.keys())
        propensities = compute_propensities(y, reactions, rate_constants)

        while time < t_end:
            total_rate = float(np.sum(propensities))
            if total_rate == 0:
                break

            # Time until the next reaction
            tau = np.random.exponential(1 / total_rate)
            time += tau

            # Determine which reaction occurs
            chosen_reaction = random.choices(
                r_keys, weights=propensities.tolist(), k=1)[0]

            # Update concentrations based on the chosen reaction
            rxn = reactions[chosen_reaction]
            species: set[int] = set()
            for reactant in rxn.reactants:
                y[reactant-1] -= 1
                species.add(reactant)
            for prod in rxn.products:
                y[prod-1] += 1
                species.add(prod)

            affected_rxns = get_affected_reactions(sp_to_reaction_map, species)
            propensities = update_affected_propensities(
                propensities, rate_constants,
                y, reactions, affected_rxns)

            self.record_line(time, y, weights, recording_times)

        return y
