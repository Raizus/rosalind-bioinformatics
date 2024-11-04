
from BioInfoToolkit.RuleBasedModel.actions.actions import BNGLACtion, GenerateNetworkAction, ResetConcentrationsAction, SaveConcentrationsAction, SetConcentrationAction, SimulateAction
from BioInfoToolkit.RuleBasedModel.model.compartment import Compartment
from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Observable import Observable
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule
from BioInfoToolkit.RuleBasedModel.model.Species import Species
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import parse_comment

def parse_block_lines(current_line: str,
                      current_block: str | None,
                      model: Model):
    molecule_types = model.molecule_types_block.items

    if current_block == 'parameters':
        parameter = Parameter.from_declaration(current_line)
        model.parameters_block.add_parameter(parameter)

    elif current_block == 'compartments':
        compartment = Compartment.from_declaration(current_line)
        model.compartments_block.add_compartment(compartment)

    elif current_block in ('species', 'seed species'):
        species = Species.from_declaration(
            current_line, molecule_types)
        model.species_block.add_species(species)

    elif current_block == 'reaction rules':
        rule = ReactionRule.from_declaration(
            current_line, molecule_types)
        model.reaction_rules_block.add_rule(rule)

    elif current_block == 'observables':
        observable = Observable.from_declaration(
            current_line, molecule_types)
        model.observables_block.add_observable(observable)


def parse_actions(current_line: str):
    action: BNGLACtion | None = None
    if current_line.startswith('generate_network'):
        action = GenerateNetworkAction.from_declaration(
            current_line)

    elif current_line.startswith('simulate'):
        action = SimulateAction.from_declaration(current_line)

    elif current_line == 'saveConcentrations();':
        action = SaveConcentrationsAction()

    elif current_line == 'resetConcentrations();':
        action = ResetConcentrationsAction()

    elif current_line.startswith('setConcentration'):
        action = SetConcentrationAction.from_declaration(current_line)

    return action


def load_bngl(file_path: str):
    model = Model()
    actions: list[BNGLACtion] = []

    current_block: str | None = None
    molecule_types_parsed: bool = False
    parsed_model: bool = False
    parsing_model: bool = False

    current_line = ""
    with open(file_path, 'r', encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            # Skip empty lines and comments
            if not line or parse_comment(line) is not None:
                continue

            # Handle line continuation (if line ends with a backslash '\')
            if line.endswith('\\'):
                current_line += line[:-1].strip() + " "
                continue
            current_line += line

            # Skip empty lines
            if not current_line:
                continue

            # Detect the beginning of a block
            if line.startswith('begin'):
                block_name = line.split(' ', 1)[1].strip()
                if block_name == 'model':
                    parsing_model = True
                current_block = block_name
                current_line = ""
                continue

            # Detect the end of a block
            if line.startswith('end'):
                block_name = line.split(' ', 1)[1].strip()
                if block_name == 'model' and parsing_model:
                    parsing_model = False
                    parsed_model = True
                current_block = None
                current_line = ""
                continue

            # Ensure molecule types block appears before other blocks
            if (current_block in ['species', 'reaction rules', 'observables']
                and not molecule_types_parsed):
                msg = ("Molecule types block must be parsed before species, "
                       + "reaction rules, or observables blocks.")
                raise ValueError(msg)

            # Parse lines depending on the current block
            if current_block == 'molecule types':
                molecule_type = MoleculeType.from_declaration(current_line)
                model.molecule_types_block.add_molecule_type(molecule_type)
                molecule_types_parsed = True  # Mark molecule types as parsed
                current_line = ""
                continue

            parse_block_lines(current_line, current_block, model)

            # parse actions
            if current_block is None:
                action = parse_actions(current_line)
                if action:
                    actions.append(action)

            # Reset current_line after processing
            current_line = ""

    model.filename = file_path
    return model, actions
