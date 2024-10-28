from collections import defaultdict
from functools import reduce
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup


import numpy as np
import numpy.typing as npt

from typing import OrderedDict

from BioInfoToolkit.RuleBasedModel.network.reaction_generation import Reaction
from BioInfoToolkit.RuleBasedModel.utils.utls import write_to_csv


def get_groups_weight_matrix(groups: OrderedDict[int, 'ObservablesGroup'], n: int):
    """
    Returns a n by m matrix (numpy array) where m is the length of groups, and 
    each column of the matrix is given by group.get_weights_array(n).

    Parameters:
    - groups: OrderedDict[int, ObservablesGroup] - a dictionary of groups
    - n: int - the length of the weights array for each group

    Returns:
    - A numpy array of shape (n, m)
    """
    # Number of groups (m) is the length of the dictionary
    m = len(groups)

    # Initialize an empty n x m matrix
    weight_matrix = np.empty((n, m), dtype=np.float64)

    # Iterate over groups and fill the matrix
    for idx, group in enumerate(groups.values()):
        weight_matrix[:, idx] = group.get_weights_array(n)

    return weight_matrix


def compute_reaction_rates(
    rate_constants: OrderedDict[int, float],
    concentrations: OrderedDict[int, int],
    reactions: OrderedDict[int, Reaction]
):
    rates: OrderedDict[int, float] = OrderedDict()
    # compute the rates for each reaction
    # the rate is equal to the product of the concentrations of the left
    # side reagents and the reaction rate constant
    for r_id, rxn in reactions.items():
        rate_constant = rate_constants[r_id]
        concent_reactants: list[int] = []
        for reactant in rxn.reactants:
            if reactant in concentrations:
                concent_reactants.append(concentrations[reactant])
            else:
                raise KeyError(
                    f"Species with id '{reactant}' not in the concentrations dictionary.")
        rate = rate_constant * \
            reduce(lambda x, y: x*y, concent_reactants, 1.0)
        rates[r_id] = rate

    return rates


def create_cdat(cdat_filename: str, n_species: int):
    """Creates a new cdat file and writes the header

    Args:
        cdat_filename (str): _description_
        n_species (int): _description_
    """
    header = ['Time'] + [f"S{i+1}" for i in range(n_species)]
    write_to_csv(cdat_filename, 'w', [header])


def create_gdat(gdat_filename: str, groups: OrderedDict[int, ObservablesGroup]):
    header = ['Time'] + [group.name for group in groups.values()]
    write_to_csv(gdat_filename, 'w', [header])


def write_data_row(filename: str, time: float, y: npt.NDArray[np.float_]):
    row = [time] + list(y)
    write_to_csv(filename, 'a', [row])


def compute_propensities(
    concentrations: npt.NDArray[np.float_],
    reactions: OrderedDict[int, Reaction],
    rate_constants: OrderedDict[int, float],
):
    propensities = np.zeros(len(reactions))

    # Calculate the propensity for each reaction
    for reaction_id, reaction in reactions.items():
        rate_constant = rate_constants[reaction_id]
        propensity = rate_constant

        # Multiply by the concentration of each reactant
        for reactant in reaction.reactants:
            if concentrations[reactant-1] <= 0:
                propensity = 0  # If any reactant's population is zero, propensity is zero
                break
            propensity *= concentrations[reactant-1]

        propensities[reaction_id] = propensity

    return propensities


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


def get_species_to_reaction_map(reactions: OrderedDict[int, Reaction]) -> defaultdict[int, set[int]]:
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


def update_affected_propensities(
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
