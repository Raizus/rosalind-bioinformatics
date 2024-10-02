
from itertools import product
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern
from BioInfoToolkit.RuleBasedModel.model.Species import Species, species_match_gen
from BioInfoToolkit.RuleBasedModel.network.blocks import GroupsBlock, ParametersBlock, ReactionsBlock, SpeciesBlock
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction


class ReactionNetwork:
    parameters_block: ParametersBlock
    species_block: SpeciesBlock
    reactions_block: ReactionsBlock
    groups_block: GroupsBlock

    graph: nx.DiGraph
    model: Model | None

    def __init__(self, model: Model | None = None) -> None:
        self.parameters_block = ParametersBlock()
        self.species_block = SpeciesBlock()
        self.reactions_block = ReactionsBlock()
        self.groups_block = GroupsBlock()
        self.model = model

    def generate_parameters(self):
        params_block = self.parameters_block
        if self.model is None:
            raise ValueError("Model must be defined in order to generate parameters.")
        model_params = self.model.parameters_block.items
        params_copy = model_params.copy()

        while len(params_copy):
            for name, param in params_copy.items():
                if name in params_block.items:
                    continue
                param2 = param.copy()
                success = param2.eval()
                if success:
                    params_block.add_parameter(param2)
                # value = eval(declaration, {"__builtins__": None}, dict_aux)

            # remove the successfuly evaluated parameters added to the parameters block
            # from params_copy. If no param is removed stop, and there are parameters 
            # which were not evaluated successfully.
            stop = True
            for name in params_block.items:
                if name in params_copy:
                    params_copy.pop(name)
                    stop = False

            if stop:
                break

    def generate_seed_species(self):
        """Generates the initial seed species, that are used to build the
        reaction network

        Raises:
            ValueError: Raised if one of the species is well defined.
        """
        if self.model is None:
            raise ValueError(
                "Model must be defined in order to generate seed species.")

        species_dict = self.model.species_block.items
        mol_types = dict(self.model.molecule_types_block.items)

        for _, specie in species_dict.items():
            if not specie.validate(mol_types):
                raise ValueError(f"Specie {specie} is not correctly defined.")
            self.species_block.add_species(specie)

    def populate_groups(self):
        """After the reactions and species have been generated,
        call this function to populate the observables groups
        """
        if self.model is None:
            raise ValueError(
                "Model must be defined in order to populate groups.")

        groups_block = self.groups_block
        obs_dict = self.model.observables_block.items
        species_dict = self.species_block.items

        for name, observable in obs_dict.items():
            group_list: list[tuple[int,int]] = []
            for sp_id, specie in species_dict.items():
                counts = observable.match_species(specie.pattern)
                if counts > 0:
                    group_list.append((sp_id, counts))

            group = ObservablesGroup(name, group_list)
            groups_block.add_group(group)

    def reaction_gen_iter(self, n_iter: int):
        reactions_block = self.reactions_block
        reactions_dict = reactions_block.items

        if self.model is None:
            raise ValueError(
                "Model must be defined in order to generate reactions.")

        species_block = self.species_block
        species_dict = species_block.items

        ni_rxs = len(reactions_dict)
        ni_species = len(species_dict)

        reaction_rules = self.model.reaction_rules_block.items

        # TODO: add memoization to make it efficient
        for rule_id, rule in reaction_rules.items():
            name = rule.name
            reactants = rule.reactants
            products = rule.products

            reactants_gens = [species_match_gen(patt, species_dict) for patt in reactants]
            # TODO: Note that, to apply a transformation to a species matching the pattern
            # we must know the node matching map, so that we can apply the transformations
            # to the correct reaction center
            # the reaction center in the reaction rule patterns may not match the center in the
            # species patterns

            for react_sp_ids in product(*reactants_gens):
                react_sp_patts = [species_dict[sp_id].pattern for sp_id in react_sp_ids]

                # apply rule to reactants
                prod_sp_patts: list[Pattern] = []

                # if generated products are new, then add them to species and 
                # add reaction to the reaction block
                prod_sp_ids: list[int] = []
                for prod_sp in prod_sp_patts:
                    specie = Species(prod_sp, "0")
                    sp_id = species_block.add_species(specie)
                    prod_sp_ids.append(sp_id)

                    # if none of the species is new, we can ignore the reaction
                    comment = f"{rule.name}"
                    rxn = Reaction(list(react_sp_ids), prod_sp_ids, rule_id, comment)
                    # create reaction

    def generate_network(self):
        print("Generate seed species")
        self.generate_seed_species()

        reactions_block = self.reactions_block
        reactions_dict = reactions_block.items

        species_block = self.species_block
        species_dict = species_block.items

        ni_rxs = len(reactions_dict)
        ni_species = len(species_dict)

        n_iter: int = 0

        while True:
            self.reaction_gen_iter(n_iter)

            ni2_rxs = len(reactions_dict)
            ni2_species = len(species_dict)

            # no new reactions or species this iteration
            if ni_rxs == ni2_rxs and ni_species == ni2_species:
                break

            n_iter += 1

