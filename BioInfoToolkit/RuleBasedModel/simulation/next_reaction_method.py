from collections import defaultdict
from typing import OrderedDict

import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import compute_propensities, \
    create_cdat, create_gdat, get_groups_weight_matrix, write_data_row
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict


def get_reactant_to_reaction_map(reactions: OrderedDict[int, Reaction]) -> defaultdict[int, set[int]]:
    # maps species_id to set of r_ids
    map_: defaultdict[int, set[int]] = defaultdict(set)

    for r_id, reaction in reactions.items():
        for reactant in reaction.reactants:
            map_[reactant].add(r_id)
    return map_


def get_affected_reactions(sp_to_reaction_map: dict[int, set[int]], species: set[int]):
    affected: set[int] = set()

    for specie in species:
        affected.update(sp_to_reaction_map[specie])

    return list(affected)


def update_affected_rates(
    propensities: npt.NDArray[np.float_],
    rate_constants: OrderedDict[int, float],
    concentrations: npt.NDArray[np.float_],
    reactions: OrderedDict[int, Reaction],
    affected: list[int]
):
    for r_id in affected:
        reaction = reactions[r_id]
        rate = rate_constants[r_id]
        for reactant in reaction.reactants:
            conc = concentrations[reactant-1]
            if conc <= 0:
                rate = 0
                break
            rate *= conc

        propensities[r_id] = rate
    return propensities


class NextReactionMethod:
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

    def simulate(
        self,
        concentrations: OrderedDict[int, int],
        reactions: OrderedDict[int, Reaction],
        rate_constants: OrderedDict[int, float],
        groups: OrderedDict[int, ObservablesGroup]
    ):

        num_species = len(concentrations)

        t_start = self.sim_params['t_start']
        t_end = self.sim_params['t_end']
        n_steps = self.sim_params['n_steps']

        # If n_steps is provided, generate time points at which to record data
        if n_steps is not None:
            recording_times = np.linspace(t_start, t_end, n_steps)
            next_recording_idx = 1  # to track the next recording time index

        time = t_start
        times: list[float] = [t_start]

        sp_to_reaction_map = get_reactant_to_reaction_map(reactions)

        y = np.array(list(concentrations.values()), dtype=np.float64)
        weights = get_groups_weight_matrix(groups, num_species)

        # write cdat file
        create_cdat(self.cdat_filename, y.size)
        write_data_row(self.cdat_filename, time, y)

        # write to gdat file
        groups_conc = y @ weights
        create_gdat(self.gdat_filename, groups)
        write_data_row(self.gdat_filename, time, groups_conc)

        propensities = compute_propensities(y, reactions, rate_constants)
        taus = -np.log(np.random.rand(propensities.size)) / propensities

        while time < t_end:
            # Find the reaction with the smallest tau (next to happen)
            next_reaction_index = int(np.argmin(taus))
            tau_min = taus[next_reaction_index]

            # Advance time by the smallest tau
            time += tau_min
            times.append(time)

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
            new_propensities = update_affected_rates(
                new_propensities, rate_constants,
                y, reactions, affected_rxns)

            # Update taus: only those reactions affected by the change need updates
            for r_id, reaction in reactions.items():
                new_rate = new_propensities[r_id]
                old_rate = propensities[r_id]
                if new_rate == 0:
                    taus[r_id] = np.inf
                elif old_rate == 0 or r_id == next_reaction_index:
                    taus[r_id] = -np.log(np.random.rand()) / new_rate
                else:
                    ratio = old_rate / new_rate
                    taus[r_id] = ratio * (taus[r_id] - tau_min)

            propensities = new_propensities

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
