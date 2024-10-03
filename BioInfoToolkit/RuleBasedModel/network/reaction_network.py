
from collections import defaultdict
import networkx as nx
import graphviz

from BioInfoToolkit.RuleBasedModel.model.Model import InvalidModelBlockError, Model
from BioInfoToolkit.RuleBasedModel.network.blocks import GroupsBlock, ParametersBlock
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import ReactionGenerator, build_rules_dict
from BioInfoToolkit.RuleBasedModel.network.reaction_block import ReactionsBlock
from BioInfoToolkit.RuleBasedModel.network.species_block import SpeciesBlock


class NetworkConstructionError(Exception):
    pass


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

    def build_network(self, model: Model):
        self.model = model

        try:
            model.validate()
        except InvalidModelBlockError as exc:
            msg = "Could not build the reaction network. Model has invalid block(s)."
            raise NetworkConstructionError(msg) from exc

        # generate parameters
        self.generate_parameters(model)

        # generate seed species
        self.generate_seed_species(model)

        # generate reactions
        self.generate_reactions(model)

        # generate groups
        self.populate_groups(model)
        a = 0

    def generate_parameters(self, model: Model):
        params_block = self.parameters_block
        model_params = model.parameters_block.items
        evaluated_params = model.parameters_block.evaluated_params

        for name in evaluated_params.keys():
            param = model_params[name]
            params_block.add_parameter(param)

    def generate_seed_species(self, model: Model):
        """Generates the initial seed species, that are used to build the
        reaction network

        Raises:
            ValueError: Raised if one of the species is well defined.
        """

        species_dict = model.species_block.items
        mol_types = dict(model.molecule_types_block.items)

        for _, specie in species_dict.items():
            if not specie.validate(mol_types):
                raise ValueError(f"Specie {specie} is not correctly defined.")
            self.species_block.add_specie(specie)

    def populate_groups(self, model: Model):
        """After the reactions and species have been generated,
        call this function to populate the observables groups
        """
        groups_block = self.groups_block
        obs_dict = model.observables_block.items
        species_dict = self.species_block.items

        for name, observable in obs_dict.items():
            group_list: list[tuple[int,int]] = []
            for sp_id, specie in species_dict.items():
                counts = observable.match_species(specie.pattern)
                if counts > 0:
                    group_list.append((sp_id, counts))

            group = ObservablesGroup(name, group_list)
            groups_block.add_group(group)

    def generate_reactions(self, model: Model):

        reactions_block = self.reactions_block
        reactions_dict = reactions_block.items

        species_block = self.species_block
        species_dict = species_block.items

        reaction_rules = build_rules_dict(model.reaction_rules_block.items)
        reaction_gen = ReactionGenerator(reaction_rules)

        n_rxs_prev = len(reactions_dict)
        n_species_prev = len(species_dict)

        n_iter: int = 0
        max_iter: int = 100

        msg = f"Iteration {n_iter}: \t {n_species_prev} species \t {n_rxs_prev} rxns"
        print(msg)

        while n_iter < max_iter:
            for reaction in reaction_gen.generate(species_block):
                self.reactions_block.add_reaction(reaction)

            n_rxs = len(reactions_dict)
            n_species = len(species_dict)

            # no new reactions or species this iteration, stop
            if n_rxs_prev == n_rxs and n_species_prev == n_species:
                break

            n_rxs_prev = n_rxs
            n_species_prev = n_species
            n_iter += 1

            msg = f"Iteration {n_iter}: \t {n_species_prev} species \t {n_rxs_prev} rxns"
            print(msg)

    def as_string(self) -> str:
        out = ''
        parameters_block_str = self.parameters_block.gen_string()
        species_block_str = self.species_block.gen_string()
        reactions_block_str = self.reactions_block.gen_string()
        groups_block_str = self.groups_block.gen_string()

        strings = [parameters_block_str,
                   species_block_str,
                   reactions_block_str,
                   groups_block_str]

        for string in strings:
            if string:
                out += string

        return out

    def draw_graph(self, filename='network_graph'):
        dot = graphviz.Digraph(comment='Reaction Network Graph', format='png')

        # maps r_id to set[(r_id, sp_id)]
        # it could be a multi digraph
        adj_dict: defaultdict[int, set[tuple[int, int]]] = defaultdict(set)
        input_species_reaction_dict: defaultdict[int, set[int]] = defaultdict(set)
        output_species_reaction_dict: defaultdict[int, set[int]] = defaultdict(set)

        for r_id, reaction in self.reactions_block.items.items():
            label = str(r_id)
            dot.node(str(r_id), label=label)
            reactants = reaction.reactants
            products = reaction.products

            for react in reactants:
                input_species_reaction_dict[react].add(r_id)
            for prods in products:
                output_species_reaction_dict[prods].add(r_id)
        
        # build adj_dict
        for sp, reactions in output_species_reaction_dict.items():
            for r_id in reactions:
                out_reactions = input_species_reaction_dict.get(sp, set())
                for r2_id in out_reactions:
                    adj_dict[r_id].add((r2_id, sp))

        for r_in, rs_out in adj_dict.items():
            for r2_id, sp in rs_out:
                # species = self.species_block.items[sp]
                label = str(sp)
                dot.edge(str(r_in), str(r2_id), label=label)

        dot.render(filename)

    def draw_hypergraph(self, filename='network_hypergraph'):
        dot = graphviz.Digraph(comment='Reaction Network Hypergraph', format='png')

        # each reaction is a hyper edge conneting species nodes

        for sp_id, species in self.species_block.items.items():
            label = str(sp_id)
            dot.node(str(sp_id), label=label)

        # build hyperedges:
        for r_id, reaction in self.reactions_block.items.items():
            node_id = f"he{r_id}"
            dot.node(node_id, shape='point', width='0.08', xlabel=f"R_{r_id}")

            for react in reaction.reactants:
                dot.edge(str(react), node_id)
            for prod in reaction.products:
                dot.edge(node_id, str(prod))

        dot.render(filename)
