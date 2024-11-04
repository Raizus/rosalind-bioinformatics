from typing import OrderedDict

from BioInfoToolkit.RuleBasedModel.model.compartment import Compartments
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern
from BioInfoToolkit.RuleBasedModel.model.Species import Species, find_species_match
from BioInfoToolkit.RuleBasedModel.network.blocks import NetworkBlock
from BioInfoToolkit.RuleBasedModel.utils.utls import eval_expr, format_data_into_lines


class SpeciesBlock(NetworkBlock):
    name = "species"
    items: OrderedDict[int, Species]
    id_count: int

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()
        self.id_count = 1

    def get_specie(self, sp_id: int) -> Species | None:
        return self.items.get(sp_id, None)

    def has_specie(self, sp_id: int) -> bool:
        specie = self.get_specie(sp_id)
        return specie is not None

    def get_specie_id(self, species_pattern: Pattern) -> int:
        """Returns the id of the species in the species dictionary that matches the 
        given species_pattern

        Args:
            target_species (Species): _description_

        Returns:
            int: species id (-1 if no match)
        """
        match_id = find_species_match(species_pattern, self.items)
        return match_id

    def add_specie(self, specie: Species):
        """Adds a new species to the block if it doesn't already exist.
        If the specie is added to the block returns its id. Else return the id of the matching species

        Args:
            species (Species): _description_

        Returns:
            _type_: _description_
        """
        species_dict = self.items
        match_id = self.get_specie_id(specie.pattern)
        if match_id in species_dict:
            return match_id

        sp_id = self.id_count
        self.items[sp_id] = specie
        self.id_count += 1
        return sp_id

    def validate_species(
        self,
        molecule_types: dict[str, MoleculeType],
        compartments: Compartments
    ) -> bool:
        for _, specie in self.items.items():
            valid = specie.validate(molecule_types, compartments)
            if not valid:
                return False
        return True

    def validate_expressions(self, variables: dict[str, int | float]) -> bool:
        for _, specie in self.items.items():
            value, _ = eval_expr(specie.expression, variables)
            if not isinstance(value, int) and not isinstance(value, float):
                return False
            specie.conc = value
        return True

    def gen_string(self):
        if len(self.items) == 0:
            return ''

        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str]] = []
        for sp_id, species in self.items.items():
            data.append((str(sp_id), str(species.pattern)))

        data_lines = format_data_into_lines(data)
        for i, (line, specie) in enumerate(zip(data_lines, self.items.values())):
            line += f" {specie.expression}"
            data_lines[i] = line
        lines.extend(data_lines)
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)
