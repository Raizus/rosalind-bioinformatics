

from BioInfoToolkit.RuleBasedModel.model.Compartment import Compartment
from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Observable import Observable
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule
from BioInfoToolkit.RuleBasedModel.model.Species import Species
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import parse_comment


def load_model(file_path: str):
    model = Model()

    current_block = None
    molecule_types_parsed = False
    molecule_types = model.molecule_types_block.items

    with open(file_path, 'r', encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Detect the beginning of a block
            if line.startswith('begin'):
                current_block = line.split(' ', 1)[1].strip()
                continue

            # Detect the end of a block
            if line.startswith('end'):
                current_block = None
                continue

            # ignore comments
            if parse_comment(line):
                continue

            # Ensure molecule types block appears before other blocks
            if (current_block in ['species', 'reaction rules', 'observables']
                and not molecule_types_parsed):
                msg = ("Molecule types block must be parsed before species, " +
                       "reaction rules, or observables blocks.")
                raise ValueError(msg)

            # Parse lines depending on the current block
            if current_block == 'molecule types':
                molecule_type = MoleculeType.from_declaration(line)
                model.molecule_types_block.add_molecule_type(molecule_type)
                molecule_types_parsed = True  # Mark molecule types as parsed

            elif current_block == 'parameters':
                parameter = Parameter.from_declaration(line)
                model.parameters_block.add_parameter(parameter)

            elif current_block == 'compartments':
                compartment = Compartment.from_declaration(line)
                model.compartments_block.add_compartment(compartment)

            elif current_block == 'species':
                species = Species.from_declaration(line, molecule_types)
                model.species_block.add_species(species)

            elif current_block == 'reaction rules':
                rule = ReactionRule.from_declaration(line, molecule_types)
                model.reaction_rules_block.add_rule(rule)

            elif current_block == 'observables':
                observable = Observable.from_declaration(line, molecule_types)
                model.observables_block.add_observable(observable)

    return model
