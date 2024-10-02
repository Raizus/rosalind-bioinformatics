

from typing import OrderedDict
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.model.Species import Species, find_species_match
from BioInfoToolkit.RuleBasedModel.utils.utls import format_data_into_lines
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction


class NetworkBlock:
    name: str

class ParametersBlock(NetworkBlock):
    name = "parameters"
    items: OrderedDict[str, Parameter]

    def __init__(self) -> None:
        self.items = OrderedDict()

    def add_parameter(self, parameter: Parameter):
        self.items[parameter.name] = parameter

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str, str]] = []
        for i, (name, param) in enumerate(self.items.items()):
            comment = f"# {param.comment}" if param.comment else ''
            data_line = (str(i), name, param.expression, comment)
            data.append(data_line)

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class SpeciesBlock(NetworkBlock):
    name = "species"
    items: OrderedDict[int, Species]
    id_count: int

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()
        self.id_count = 1

    def add_species(self, species: Species):
        """Adds a new species to the block if it doesn't already exist.
        If the specie is added to the block returns it's id. Else return the id of the matching species

        Args:
            species (Species): _description_

        Returns:
            _type_: _description_
        """
        species_dict = self.items
        match_id = find_species_match(species.pattern, species_dict)
        if match_id in species_dict:
            return match_id

        sp_id = self.id_count
        self.items[sp_id] = species
        self.id_count += 1
        return sp_id

    def validate_species(self, molecule_types: dict[str, MoleculeType]) -> bool:
        for _, specie in self.items.items():
            valid = specie.validate(molecule_types)
            if not valid:
                return False
        return True

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str]] = []
        for sp_id, species in self.items.items():
            data.append((str(sp_id), str(species.pattern), species.expression))

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class ReactionsBlock(NetworkBlock):
    name = "reactions"
    items: OrderedDict[int, Reaction]
    count_id: int

    def __init__(self) -> None:
        self.count_id = 1
        self.items = OrderedDict()

    def add_reaction(self, reaction: Reaction):
        r_id = self.count_id
        self.items[r_id] = reaction
        self.count_id += 1

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str, str, str]] = []
        for r_id, reaction in self.items.items():
            comment = f"# {reaction.comment}" if reaction.comment else ''
            data_line = (str(r_id),
                         str(reaction.reactants),
                         str(reaction.products),
                         str(reaction.rate_expression),
                         comment)
            data.append(data_line)

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class GroupsBlock(NetworkBlock):
    name = "groups"
    items: OrderedDict[int, ObservablesGroup]
    count_id: int

    def __init__(self) -> None:
        self.items = OrderedDict()
        self.count_id = 1

    def add_group(self, group: ObservablesGroup):
        self.items[self.count_id] = group
        self.count_id += 1

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str]] = []
        for g_id, group in self.items.items():
            group_str = ','.join(f'{w}*{sp_id}' if w != 1 else f'{sp_id}' for sp_id,
                                 w in group.weighted_species)
            data.append((str(g_id), str(group.name), group_str))

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)
