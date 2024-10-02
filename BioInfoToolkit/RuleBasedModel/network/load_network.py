

# from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
# from BioInfoToolkit.RuleBasedModel.model.Species import Species
# from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
# from BioInfoToolkit.RuleBasedModel.network.parsers import parse_parameters, parse_reaction, parse_seed_species
# from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
# from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork

# def parse_net_file(file):
#     # parse file


# def load_network(file_path: str):
#     network = ReactionNetwork()

#     # load .net file in path

#     # if parameter block
#     parsed_parameter = parse_parameters(line)
#     name = parsed_parameter['name']
#     expression = parsed_parameter['expression']
#     parameter = Parameter(name, expression)
#     network.parameters_block.add_parameter(parameter)

#     # if species block
#     parsed_species = parse_seed_species(line)
#     species = Species.
#     network.species_block.add_species()

#     # if reactions block
#     reaction = Reaction.from_declaration(line)
#     network.reactions_block.add_reaction(reaction)

#     # if groups block
#     group = ObservablesGroup.from_declaration(line)
#     network.groups_block.add_group(group)
