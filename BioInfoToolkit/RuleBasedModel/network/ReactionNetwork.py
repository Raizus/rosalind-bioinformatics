
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule
from BioInfoToolkit.RuleBasedModel.model.Species import Species
from BioInfoToolkit.RuleBasedModel.network.Reaction import Reaction


class ReactionNetwork:
    graph: nx.DiGraph
    reaction_rules: list[ReactionRule]
    species: dict[int, Species]
    reactions: list[Reaction]

    def __init__(self, reaction_rules: list[ReactionRule], seed_species: list[Pattern]) -> None:
        # Each node represents a reaction (maybe reaction rule in general?)
        # each outgoing edge corresponds to a reaction product
        # each incoming edge corresponds to a reaction reactant
        # should each node be a single order transformation or correspond to a full rule?
        # reactants and products in the rules must be well defined species.
        species: dict[int, Species] = {}
        self.species = species
        for i, pattern in enumerate(seed_species, start=1):
            specie = Species(pattern, i)
            # TODO: make sure specie is not duplicate
            species[i] = specie

        # TODO: remove reactions with zero rates
        self.reaction_rules = reaction_rules
        self.reactions = []

    def reaction_gen_iter(self):
        rxs_l = len(self.reactions)
        species_l = len(self.species)

        for rule in self.reaction_rules:
            # find all matches between reactant patterns and species
            

            # apply rule to all possible combinations to generate reactions

            # we need only to consider reactions that involve new species that were generated in the last iteration

            # add reaction if not already in the list
            pass
        
        # if no new species / reactions then we can finish iterating