
import os
from collections import defaultdict
from typing import OrderedDict
import networkx as nx
import graphviz
import numpy as np

from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.model.ModelBlocks import InvalidModelBlockError
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.model.Species import Species
from BioInfoToolkit.RuleBasedModel.network.blocks import GroupsBlock, ParametersBlock
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction_generation import Reaction, ReactionGenerator, \
    build_rules_dict
from BioInfoToolkit.RuleBasedModel.network.reaction_block import ReactionsBlock
from BioInfoToolkit.RuleBasedModel.network.species_block import SpeciesBlock
from BioInfoToolkit.RuleBasedModel.simulation.gillespie import GillespieSimulator
from BioInfoToolkit.RuleBasedModel.simulation.next_reaction_method import NextReactionMethod
from BioInfoToolkit.RuleBasedModel.simulation.ode_sim import ODESimulator
from BioInfoToolkit.RuleBasedModel.simulation.progressive_leaping import ProgressiveLeapingSimulator
from BioInfoToolkit.RuleBasedModel.simulation.tau_leaping import TauLeapingSimulator
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import GenerateNetworkDict, SimulateDict
from BioInfoToolkit.RuleBasedModel.utils.network_parsers import parse_parameters, parse_seed_species
from BioInfoToolkit.RuleBasedModel.utils.utls import eval_expr, compose_path, decompose_path


class NetworkConstructionError(Exception):
    pass


class ReactionNetwork:
    parameters_block: ParametersBlock
    species_block: SpeciesBlock
    reactions_block: ReactionsBlock
    groups_block: GroupsBlock
    concentrations: OrderedDict[int, int]

    graph: nx.DiGraph
    model: Model | None
    net_filename: str
    cdat_filename: str
    gdat_filename: str

    def __init__(self) -> None:
        self.parameters_block = ParametersBlock()
        self.species_block = SpeciesBlock()
        self.reactions_block = ReactionsBlock()
        self.groups_block = GroupsBlock()
        self.model = None
        self.concentrations = OrderedDict()

        self.net_filename = 'model.net'
        self.cdat_filename = 'model.cdat'
        self.gdat_filename = 'model.gdat'

    def load_network(self, net_path: str):
        network = load_network(net_path)
        self.parameters_block = network.parameters_block
        self.reactions_block = network.reactions_block
        self.species_block = network.species_block
        self.groups_block = network.groups_block
        self.set_output_filepaths(net_path)

    def generate_network(self, model: Model,
                         params: GenerateNetworkDict):
        overwrite = params['overwrite']
        path, name, _ = decompose_path(model.filename)
        net_path = compose_path(path, name, ".net")
        # if overwrite is False, and a .net file already exists, load that file instead
        if not overwrite and os.path.isfile(net_path):
            self.load_network(net_path)
            return

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
        max_iter = params['max_iter']
        max_stoich = params['max_stoich']
        self.generate_reactions(model, max_iter, max_stoich)

        # generate groups
        self.populate_groups(model)

        self.set_output_filepaths(net_path)

    def set_output_filepaths(self, net_filename: str):
        """Given the filepath of the net file, sets the net_filename, cdat_filename 
        and gdat_filename, with the same path and name as the net file

        Args:
            net_filename (str): _description_
        """
        path, name, ext = decompose_path(net_filename)
        net_path = compose_path(path, name, ext)
        self.net_filename = net_path
        cdat_path = compose_path(path, name, ".cdat")
        self.cdat_filename = cdat_path
        gdat_path = compose_path(path, name, ".gdat")
        self.gdat_filename = gdat_path

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
        compartments = model.compartments_block.as_compartments()

        for _, specie in species_dict.items():
            if not specie.validate(mol_types, compartments):
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
            # maps species id to weights
            group_dict: defaultdict[int, int] = defaultdict(int)
            for sp_id, specie in species_dict.items():
                counts = observable.match_species(specie.pattern)
                if counts == 0:
                    continue
                if observable.type == 'Molecules':
                    group_dict[sp_id] += counts
                elif observable.type == 'Species':
                    group_dict[sp_id] = counts

            group = ObservablesGroup(name, group_dict)
            groups_block.add_group(group)

    def generate_reactions(self, model: Model,
                           max_iter: int | None = None,
                           max_stoich: dict[str, int] | None = None):

        reactions_block = self.reactions_block
        reactions_dict = reactions_block.items

        species_block = self.species_block
        species_dict = species_block.items

        reaction_rules = build_rules_dict(model.reaction_rules_block.items)
        # substitute reaction expression if they're not variables
        r_law_id = 1
        for _, reaction in reaction_rules.items():
            expr = reaction.forward_rate
            if expr in self.parameters_block.items:
                continue

            while True:
                var_name = f"_rateLaw{r_law_id}"
                if var_name not in self.parameters_block.items:
                    break
                r_law_id += 1
            new_param = Parameter(var_name, expr)
            self.parameters_block.add_parameter(new_param)
            reaction.forward_rate = var_name
            r_law_id += 1

        reaction_gen = ReactionGenerator(reaction_rules)

        n_rxs_prev = len(reactions_dict)
        n_species_prev = len(species_dict)

        n_iter: int = 0
        if max_stoich is None:
            max_stoich = {}

        msg = f"Iteration {n_iter}:\t{n_species_prev} species\t{n_rxs_prev} rxns"
        print(msg)

        while True:
            for reaction in reaction_gen.generate(species_block, max_stoich):
                self.reactions_block.add_reaction(reaction)

            n_rxs = len(reactions_dict)
            n_species = len(species_dict)

            # no new reactions or species this iteration, stop
            if n_rxs_prev == n_rxs and n_species_prev == n_species:
                break

            n_rxs_prev = n_rxs
            n_species_prev = n_species
            n_iter += 1

            msg = f"Iteration {n_iter}:\t{n_species_prev} species\t{n_rxs_prev} rxns"
            print(msg)

            if max_iter and n_iter >= max_iter:
                break

        msg = f"TOTAL {n_iter}:\t{n_species_prev} species\t{n_rxs_prev} rxns"
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

        for sp_id, _ in self.species_block.items.items():
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

    def get_rate_constants(self) -> OrderedDict[int, float]:
        evaluated_params = self.parameters_block.evaluated_params

        rate_constants: OrderedDict[int, float] = OrderedDict()
        for r_id, rxn in self.reactions_block.items.items():
            rate_expr = rxn.rate_expression
            rate_val, _ = eval_expr(rate_expr, evaluated_params)
            if not isinstance(rate_val, int) and not isinstance(rate_val, float):
                raise TypeError(f"rate_val '{rate_val}' must be int or float.")
            rate_constants[r_id] = float(rate_val)

        return rate_constants

    def simulate(self, params: SimulateDict):
        method = params['method']

        self.parameters_block.evaluate_parameters()

        # maps reaction id's to reaction rate constants
        rate_constants = self.get_rate_constants()
        groups = self.groups_block.items
        reactions = self.reactions_block.items
        concentrations = self.initialise_concentrations()

        if method == 'ssa':
            yf = self.gillespie_simulation(params)
        elif method == 'ode':

            # evaluate species expressions and initialize concentrations
            y = np.array(list(concentrations.values()), dtype=np.float64)

            simulator = ODESimulator(params, self.cdat_filename, self.gdat_filename)
            yf = simulator.simulate(y, reactions, rate_constants, groups)
        elif method == 'tau-leap':
            simulator = TauLeapingSimulator(params, self.cdat_filename, self.gdat_filename)
            yf = simulator.solve(concentrations, reactions, rate_constants, groups)
        elif method == 'nrm':
            simulator = NextReactionMethod(params, self.cdat_filename, self.gdat_filename)
            yf = simulator.simulate(concentrations, reactions, rate_constants, groups)
        elif method == 'pla':
            simulator = ProgressiveLeapingSimulator(
                params, self.cdat_filename, self.gdat_filename)
            yf = simulator.simulate(
                concentrations, reactions, rate_constants, groups)
        else:
            raise ValueError(f"Simulation method '{method}' is not valid.")

        for i, conc in enumerate(yf, start=1):
            self.concentrations[i] = int(conc)

    def initialise_concentrations(self):
        # evaluate species expressions and initialize concentrations
        evaluated_params = self.parameters_block.evaluated_params
        self.species_block.validate_expressions(evaluated_params)
        concentrations: OrderedDict[int, int] = OrderedDict()
        for sp_id, specie in self.species_block.items.items():
            concentrations[sp_id] = int(specie.conc)

        self.concentrations = concentrations
        return concentrations

    def reset_concentrations(self):
        self.initialise_concentrations()

    # def set_concentration(self, specie_q: int | Pattern, value: str):
    #     if isinstance(specie_q, int):
    #         specie = self.species_block.get_specie(specie_q)
    #         if specie:
    #             specie.expression = value
    #     else:
    #         pass

    def gillespie_simulation(self, params: SimulateDict):

        # evaluate species expressions and initialize concentrations
        concentrations = self.initialise_concentrations()

        # maps reaction id's to reaction rate constants
        rate_constants = self.get_rate_constants()

        reactions = self.reactions_block.items
        groups = self.groups_block.items

        simulator = GillespieSimulator(params, self.cdat_filename, self.gdat_filename)
        y = simulator.simulate(reactions, rate_constants, concentrations, groups)
        return y

    def save_concentrations(self):
        """
        Saves the current concentrations to the species block and saves the network file
        """
        concentrations = self.concentrations
        species_dict = self.species_block.items
        for i, conc in enumerate(concentrations, start=1):
            species_dict[i].conc = conc
            species_dict[i].expression = str(conc)

        self.save_network(None, True)

    def save_network(self, fp: str | None = None, overwrite: bool | None = None):
        if fp is None:
            fp = self.net_filename

        # Check if the file exists and raise an error if overwrite is False or None
        if not overwrite and os.path.exists(fp):
            msg = (f"A file already exists at '{fp}'. To overwrite, set 'overwrite' to True "
                   + "or pass overwrite=>1 to the generate_network action.")
            print(msg)
            return
            # raise FileExistsError(msg)

        with open(fp, 'w', encoding='utf-8') as net_file:
            out = self.as_string()
            net_file.write(out)


def load_network(file_path: str) -> ReactionNetwork:
    network = ReactionNetwork()

    current_block = None
    with open(file_path, 'r', encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Detect the beginning of a block
            if line.startswith('begin'):
                current_block = line.split()[1]
                continue

            # Detect the end of a block
            if line.startswith('end'):
                current_block = None
                continue

            # Parse lines depending on the current block
            if current_block == 'parameters':
                parsed_parameter = parse_parameters(line)
                parameter = Parameter.from_dict(parsed_parameter)
                network.parameters_block.add_parameter(parameter)

            elif current_block == 'species':
                parsed_species = parse_seed_species(line)
                species = Species.from_dict(parsed_species, None)
                network.species_block.add_specie(species)

            elif current_block == 'reactions':
                reaction = Reaction.from_declaration(line)
                network.reactions_block.add_reaction(reaction)

            elif current_block == 'groups':
                group = ObservablesGroup.from_declaration(line)
                network.groups_block.add_group(group)

    # set filepaths for output files
    network.set_output_filepaths(file_path)

    return network
