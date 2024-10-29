
from collections import defaultdict
from typing import OrderedDict
from BioInfoToolkit.RuleBasedModel.model.Species import Species
from BioInfoToolkit.RuleBasedModel.network.reaction_generation import Reaction


def compare_species_dicts(
    dict1: OrderedDict[int, Species],
    dict2: OrderedDict[int, Species]
) -> dict[int, int] | None:
    # defaultdict to store possible mappings
    mapping: defaultdict[int, set[int]] = defaultdict(set)
    in_cdomain: set[int] = set()

    # Populate the mapping of species in dict1 to dict2
    for sp_id1, specie1 in dict1.items():
        for sp_id2, specie2 in dict2.items():
            if specie1.pattern == specie2.pattern:
                mapping[sp_id1].add(sp_id2)
                in_cdomain.add(sp_id2)

    # Verify that every species in dict1 maps to exactly one species in dict2
    bijective_map: dict[int, int] = {}

    # Check if the map is bijective
    used_in_cdomain = set()

    for sp_id1, candidates in mapping.items():
        # Each species in dict1 must map to exactly one species in dict2
        if len(candidates) != 1:
            return None

        sp_id2 = list(candidates)[0]

        # Ensure no two species from dict1 map to the same species in dict2
        if sp_id2 in used_in_cdomain:
            return None

        bijective_map[sp_id1] = sp_id2
        used_in_cdomain.add(sp_id2)

    # Ensure every species in dict2 is accounted for
    if len(used_in_cdomain) != len(dict2):
        return None

    return bijective_map


def compare_reactions_dicts(
    dict1: OrderedDict[int, Reaction],
    dict2: OrderedDict[int, Reaction],
    species_biject_map: dict[int, int]
) -> dict[int, int] | None:
    # defaultdict to store possible mappings
    mapping: defaultdict[int, set[int]] = defaultdict(set)
    in_cdomain: set[int] = set()

    # Populate the mapping of species in dict1 to dict2
    for r_id1, reaction1 in dict1.items():
        for r_id2, reaction2 in dict2.items():
            mapped_reactants = [species_biject_map[r] for r in reaction1.reactants]
            mapped_products = [species_biject_map[r]
                                for r in reaction1.products]
            mapped_reaction = Reaction(mapped_reactants, mapped_products,
                                       reaction1.rate_expression,
                                       reaction1.rule_id, reaction1.comment)

            if mapped_reaction == reaction2:
                mapping[r_id1].add(r_id2)
                in_cdomain.add(r_id2)

    # Verify that every species in dict1 maps to exactly one species in dict2
    bijective_map: dict[int, int] = {}

    # Check if the map is bijective
    used_in_cdomain = set()

    for r_id1, candidates in mapping.items():
        # Each species in dict1 must map to exactly one species in dict2
        if len(candidates) != 1:
            return None

        r_id2 = list(candidates)[0]

        # Ensure no two species from dict1 map to the same species in dict2
        if r_id2 in used_in_cdomain:
            return None

        bijective_map[r_id1] = r_id2
        used_in_cdomain.add(r_id2)

    # Ensure every species in dict2 is accounted for
    if len(used_in_cdomain) != len(dict2):
        return None

    return bijective_map
