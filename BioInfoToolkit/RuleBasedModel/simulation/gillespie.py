from typing import OrderedDict
import random

import numpy as np

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import compute_propensities, \
    create_cdat, create_gdat, get_groups_weight_matrix, write_data_row
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict


class GillespieSimulator:
    sim_params: SimulateDict
    cdat_filename: str = 'output.cdat'
    gdat_filename: str = 'output.gdat'

    def __init__(self,
                 sim_params: SimulateDict,
                 cdat_filename: str = 'output.cdat',
                 gdat_filename: str = 'output.gdat'
                ) -> None:
        self.sim_params = sim_params
        self.cdat_filename = cdat_filename
        self.gdat_filename = gdat_filename

    def simulate(self,
                 reactions: OrderedDict[int, Reaction],
                 rate_constants: OrderedDict[int, float],
                 concentrations: OrderedDict[int, int],
                 groups: OrderedDict[int, ObservablesGroup]):
        t_start = self.sim_params['t_start']
        t_end = self.sim_params['t_end']
        time = t_start
        times: list[float] = [t_start]
        n_steps = self.sim_params['n_steps']

        # Initialize observable groups concentrations
        y = np.array(list(concentrations.values()), dtype=np.float64)
        num_species = len(concentrations)
        weights = get_groups_weight_matrix(groups, num_species)

        # If n_steps is provided, generate time points at which to record data
        if n_steps is not None:
            recording_times = np.linspace(t_start, t_end, n_steps)
            next_recording_idx = 1  # to track the next recording time index

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

        while time < t_end:
            propensities = compute_propensities(y, reactions, rate_constants)
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
            for reactant in rxn.reactants:
                y[reactant-1] -= 1
            for prod in rxn.products:
                y[prod-1] += 1

            # If n_steps is provided, only record concentrations at predefined time points
            if n_steps is not None:
                while next_recording_idx < n_steps and time >= recording_times[next_recording_idx]:
                    times.append(recording_times[next_recording_idx])

                    groups_conc = y @ weights

                    # Record concentrations of reactants
                    write_data_row(self.cdat_filename, time, y)

                    # Record concentrations of observables
                    write_data_row(self.gdat_filename, time, groups_conc)

                    print(f"\tt = {time}")

                    next_recording_idx += 1
            else:
                # Record time and concentrations for all reactions
                times.append(time)
                groups_conc = y @ weights

                # Record concentrations of reactants
                write_data_row(self.cdat_filename, time, y)

                # Record concentrations of observables
                write_data_row(self.gdat_filename, time, groups_conc)

        return times, groups_conc
