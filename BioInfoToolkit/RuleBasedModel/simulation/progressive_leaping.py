from typing import OrderedDict

import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import calculate_tau_full, compute_propensities, compute_stoichiometry_matrix, \
    create_cdat, create_gdat, get_groups_weight_matrix, write_data_row
from BioInfoToolkit.RuleBasedModel.simulation.simulator import SimulatorABC


def delta_y_poisson(
    y: npt.NDArray[np.float_],
    reactions: OrderedDict[int, Reaction],
    reaction_firings: npt.NDArray[np.int_],
    poisson_group: list[int]
):
    delta_y = np.zeros(y.shape, np.int64)

    # Update the concentrations based on how many times each reaction occurs
    for r_id, num_firings in zip(poisson_group, reaction_firings):
        reaction = reactions[r_id]

        # Apply the adjusted number of firings to the concentrations
        for reactant in reaction.reactants:
            delta_y[reactant-1] -= num_firings

        for product in reaction.products:
            delta_y[product-1] += num_firings

    return delta_y


def delta_y_langevin(
    y: npt.NDArray[np.float_],
    reactions: OrderedDict[int, Reaction],
    reaction_firings: npt.NDArray[np.float_],
    langevin_group: list[int]
):
    delta_y = np.zeros(y.shape, np.int64)

    # Update the concentrations based on how many times each reaction occurs
    for r_id, num_firings in zip(langevin_group, reaction_firings):
        reaction = reactions[r_id]

        # Apply the adjusted number of firings to the concentrations
        for reactant in reaction.reactants:
            delta_y[reactant-1] -= num_firings

        for product in reaction.products:
            delta_y[product-1] += num_firings

    return delta_y


def delta_y_ode(
    y: npt.NDArray[np.float_],
    reactions: OrderedDict[int, Reaction],
    reaction_firings: npt.NDArray[np.float_],
    ode_group: list[int]
):
    delta_y = np.zeros(y.shape, np.int64)

    # Update the concentrations based on how many times each reaction occurs
    for r_id, num_firings in zip(ode_group, reaction_firings):
        reaction = reactions[r_id]

        # Apply the adjusted number of firings to the concentrations
        for reactant in reaction.reactants:
            delta_y[reactant-1] -= num_firings

        for product in reaction.products:
            delta_y[product-1] += num_firings

    return delta_y


def classify_reactions(propensities: npt.NDArray[np.float_], tau: float):
    # devide reactions into four groups
    # if a_μ*τ <∼ 1 → Exact Stochastic(very slow)
    # if a_μ*τ > 1 but !≫ 1 → Poisson (slow)
    # If a_μ*τ ≫ 1 but √a_μ*τ  !≫ 1 → Langevin(medium)
    # If √a_μ*τ ≫ 1 → Deterministic (fast)
    es_group: list[int] = []
    poisson_group: list[int] = []
    langevin_group: list[int] = []
    deterministic_group: list[int] = []

    threshold1 = 1
    threshold2 = 100

    a_tau = propensities * tau
    for r_id, val in enumerate(a_tau):
        if val <= threshold1:   # exact stochastic
            es_group.append(r_id)
        elif threshold1 < val < threshold2: # poisson
            poisson_group.append(r_id)
        elif (val)**(1/2) < threshold2 <= val:  # langevin
            langevin_group.append(r_id)
        else:   # deterministic / ode
            # langevin_group.append(r_id)
            deterministic_group.append(r_id)

    return es_group, poisson_group, langevin_group, deterministic_group


class ProgressiveLeapingSimulator(SimulatorABC):

    def simulate(self,
              concentrations: OrderedDict[int, int],
              reactions: OrderedDict[int, Reaction],
              rate_constants: OrderedDict[int, float],
              groups: OrderedDict[int, ObservablesGroup]):

        num_species = len(concentrations)
        epsilon = 0.03

        t_start = self.sim_params['t_start']
        t_end = self.sim_params['t_end']
        n_steps = self.sim_params['n_steps']
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

        print("Starting simulation with progressive-leaping algorithm...")
        print(f"\tt = {time}")

        propensities = compute_propensities(y, reactions, rate_constants)
        taus_es = -np.log(np.random.rand(propensities.size)) / propensities

        # Perform tau-leaping simulation
        while time < t_end:
            propensities = compute_propensities(y, reactions, rate_constants)

            tau = calculate_tau_full(
                y, reactions, propensities, stoichiometry_matrix, epsilon=epsilon)

            taus_es = - \
                np.log(np.random.rand(propensities.size)) / propensities

            # tau selection loop
            while True:
                es_group, poisson_group, langevin_group, ode_group = classify_reactions(
                    propensities, tau)

                if len(es_group) > 0:
                    idx = int(np.argmin(taus_es[es_group]))
                    next_rxn_idx = es_group[idx]
                    min_tau_es = float(taus_es[next_rxn_idx])
                    if len(es_group) == len(reactions) and tau != min_tau_es:
                        tau = min_tau_es
                    elif min_tau_es < tau:
                        tau = min_tau_es
                        continue

                # compute new concentrations
                # Exact stochastic group
                dy_es = np.zeros(y.shape, np.int64)
                if len(es_group) > 0 and tau == min_tau_es:  # fire 1 reaction
                    reaction = reactions[next_rxn_idx]
                    for reactant in reaction.reactants:
                        dy_es[reactant-1] -= 1
                    for prod in reaction.products:
                        dy_es[prod-1] += 1

                # Poisson group
                # Draw the number of times each reaction occurs from a Poisson distribution
                poisson_firings = np.random.poisson(
                    propensities[poisson_group] * tau)
                dy_poisson = delta_y_poisson(
                    y, reactions, poisson_firings, poisson_group)

                # langevin group
                # Draw the number of times each reaction occurs from a Normal distribution
                a_j_t = propensities[langevin_group] * tau
                langevin_firings = a_j_t + \
                    np.sqrt(a_j_t) * np.random.rand(len(langevin_group))
                dy_langevin = delta_y_langevin(
                    y, reactions, langevin_firings, langevin_group)

                ode_firings = propensities[ode_group] * tau
                dy_ode = delta_y_ode(
                    y, reactions, ode_firings, ode_group)

                # Update the concentrations based on how many times each reaction occurs
                delta_y = np.round(dy_es + dy_poisson + dy_langevin + dy_ode)
                new_y = y + delta_y

                # check if concentrations are valid
                if not np.any(new_y < 0):
                    break

                tau = tau / 2

            # Update time
            time += tau
            y = new_y

            self.record_line(time, y, weights, recording_times)

        return y
