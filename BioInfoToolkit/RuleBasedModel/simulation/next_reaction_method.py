from collections import defaultdict
from typing import OrderedDict

import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import compute_reaction_rates, get_groups_weight_matrix
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict
from BioInfoToolkit.RuleBasedModel.utils.utls import write_to_csv


def get_reactant_to_reaction_map(reactions: OrderedDict[int, Reaction]) -> defaultdict[int, set[int]]:
    # maps species_id to set of r_ids
    map_: defaultdict[int, set[int]] = defaultdict(set)

    for r_id, reaction in reactions.items():
        for reactant in reaction.reactants:
            map_[reactant].add(r_id)
    return map_


def get_affected_reactions(sp_to_reaction_map: dict[int, set[int]], species: list[int]):
    affected: set[int] = set()

    for specie in species:
        affected.update(sp_to_reaction_map[specie])

    return list(affected)


def update_affected_rates(
    rates: OrderedDict[int, float],
    rate_constants: OrderedDict[int, float],
    concentrations: OrderedDict[int, int],
    reactions: OrderedDict[int, Reaction],
    affected: list[int]
):
    for r_id in affected:
        reaction = reactions[r_id]
        rate = rate_constants[r_id]
        for reactant in reaction.reactants:
            if concentrations[reactant] <= 0:
                rate = 0
                break
            rate *= concentrations[reactant]

        rates[r_id] = rate
    return rates


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
        time = t_start
        times: list[float] = [t_start]

        sp_to_reaction_map = get_reactant_to_reaction_map(reactions)

        n_steps = self.sim_params['n_steps']
        weights = get_groups_weight_matrix(groups, num_species)

        rates = compute_reaction_rates(
            rate_constants, concentrations, reactions)

        taus = -np.log(np.random.rand(len(rates))) \
            / np.array(rates.values(), dtype=np.float64)

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
                concentrations[reactant] -= 1
                species.add(reactant)
            for prod in reaction.products:
                concentrations[prod] += 1
                species.add(prod)

            affected_rxns = get_affected_reactions(sp_to_reaction_map, list(species))

            # Recompute the rates after the reaction affects the system
            rates = update_affected_rates(
                rates, rate_constants,
                concentrations, reactions, affected_rxns)

            # Update taus: only those reactions affected by the change need updates
            for r_id, reaction in reactions.items():
                if r_id in affected_rxns:
                    taus[r_id] = ((taus[r_id] - tau_min)
                                  if taus[r_id] > tau_min else 0)
                    if taus[r_id] == 0:
                        taus[r_id] = -np.log(np.random.rand()) / rates[r_id]
