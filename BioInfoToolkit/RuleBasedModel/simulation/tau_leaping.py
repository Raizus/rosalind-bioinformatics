from typing import OrderedDict

import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import compute_propensities, \
    create_cdat, create_gdat, get_groups_weight_matrix, write_data_row
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict


def compute_stoichiometry_matrix(reactions: OrderedDict[int, Reaction], num_species: int):
    # Stoichiometry matrix (species x reactions)
    num_reactions = len(reactions)
    stoichiometry_matrix = np.zeros((num_species, num_reactions), dtype=np.int64)

    for reaction_id, reaction in reactions.items():
        for reactant in reaction.reactants:
            stoichiometry_matrix[reactant-1, reaction_id] -= 1
        for product in reaction.products:
            stoichiometry_matrix[product-1, reaction_id] += 1

    return stoichiometry_matrix


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


def calculate_tau_full(concentrations,
                       reactions: OrderedDict[int, Reaction],
                       propensities: npt.NDArray[np.float_],
                       stoichiometry_matrix: npt.NDArray[np.int_],
                       epsilon=0.03):
    # Helper function to compute tau using the full Cao method
    num_species = len(concentrations)

    # Initialize mu_i and sigma_i^2
    mu = np.zeros(num_species)
    sigma2 = np.zeros(num_species)

    # Calculate mu_i and sigma_i^2 for each species
    for reaction_id, reaction in reactions.items():
        propensity = propensities[reaction_id]

        # Update mu_i and sigma_i^2 for each reactant and product
        for i in range(num_species):
            v_ij = stoichiometry_matrix[i, reaction_id]
            mu[i] += v_ij * propensity
            sigma2[i] += (v_ij ** 2) * propensity

    # Determine the highest-order event (g_i) for each species
    g = np.ones(num_species)  # Assuming unimolecular (1) or bimolecular (2)
    for reaction in reactions.values():
        if len(reaction.reactants) == 2:  # Bimolecular reaction
            for reactant in reaction.reactants:
                g[reactant-1] = max(g[reactant-1], 2)

    # Calculate tau
    tau = float('inf')
    for i in range(num_species):
        if mu[i] != 0:
            term1 = (max(epsilon * concentrations[i] / g[i], 1)) / abs(mu[i])
            tau = min(tau, term1)
        if sigma2[i] != 0:
            term2 = (
                max(epsilon * concentrations[i] / g[i], 1)) ** 2 / sigma2[i]
            tau = min(tau, term2)

    return tau


class TauLeapingSimulator:
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

    def solve(self,
              concentrations: OrderedDict[int, int],
              reactions: OrderedDict[int, Reaction],
              rate_constants: OrderedDict[int, float],
              groups: OrderedDict[int, ObservablesGroup]):

        num_species = len(concentrations)
        epsilon = 0.03

        y = np.array(list(concentrations.values()), dtype=np.float64)
        t_start = self.sim_params['t_start']
        time = t_start
        times: list[float] = [t_start]
        t_end = self.sim_params['t_end']
        n_steps = self.sim_params['n_steps']
        results = [y.copy()]

        weights = get_groups_weight_matrix(groups, num_species)

        # If n_steps is provided, generate time points at which to record data
        if n_steps is not None:
            recording_times = np.linspace(t_start, t_end, n_steps)
            next_recording_idx = 1  # to track the next recording time index

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

        given_tau = self.sim_params['tau']

        # Perform tau-leaping simulation
        while time < t_end:
            # Calculate propensities for each reaction
            propensities = compute_propensities(y, reactions, rate_constants)

            # Calculate tau using the full Cao method
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

            # If n_steps is provided, only record concentrations at predefined time points
            if n_steps is not None:
                while next_recording_idx < n_steps and time >= recording_times[next_recording_idx]:
                    times.append(recording_times[next_recording_idx])

                    groups_conc = y @ weights

                    # Record concentrations of reactants
                    write_data_row(self.cdat_filename, time, y)

                    # Record concentrations of observables
                    write_data_row(self.gdat_filename, time, groups_conc)

                    next_recording_idx += 1
            else:
                # Record time and concentrations for all reactions
                times.append(time)
                groups_conc = y @ weights

                # Record concentrations of reactants
                write_data_row(self.cdat_filename, time, y)

                # Record concentrations of observables
                write_data_row(self.gdat_filename, time, groups_conc)

            # Store results
            results.append(y.copy())
            print(f"\tt = {time}")

        return np.array(times), np.array(results)
