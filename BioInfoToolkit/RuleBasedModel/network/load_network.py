
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.model.Species import Species
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.utils.network_parsers import parse_parameters, parse_seed_species
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork


def load_network(file_path: str):
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
                network.species_block.add_species(species)

            elif current_block == 'reactions':
                reaction = Reaction.from_declaration(line)
                network.reactions_block.add_reaction(reaction)

            elif current_block == 'groups':
                group = ObservablesGroup.from_declaration(line)
                network.groups_block.add_group(group)

    return network
