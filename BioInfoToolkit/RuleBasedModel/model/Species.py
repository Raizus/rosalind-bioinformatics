
from itertools import product
from BioInfoToolkit.RuleBasedModel.model.Component import components_gen
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Parsers import parse_pattern
from BioInfoToolkit.RuleBasedModel.model.Pattern import Molecule, Pattern, match_patterns


class Species:
    id: int
    compartment: str|None
    pattern: Pattern

    def __init__(self, pattern: Pattern, id: int) -> None:
        # if not pattern.is_specie():
        #     raise ValueError(f"Reactant {pattern} must be a specie and not a pattern.")
        self.reactant = pattern
        self.compartment = None
        self.id = id

    @classmethod
    def from_declaration(cls,
                         declaration: str,
                         molecules: dict[str, MoleculeType],
                         specie_id: int):
        # parse string declaration
        parsed = parse_pattern(declaration)
        if not parsed:
            raise ValueError(f"Invalid reactant declaration: {declaration}")

        pattern = Pattern.from_dict(parsed, molecules)
        specie = Species(pattern, specie_id)
        return specie

    def __repr__(self) -> str:
        out = ""
        if self.compartment:
            out = f"@{self.compartment}::"
        out += str(self.reactant)
        return out

    def match_pattern(self, pattern: Pattern) -> bool:
        """Matched the specie to a given pattern

        Args:
            pattern (Pattern): _description_

        Returns:
            bool: _description_
        """
        species_pattern = self.pattern
        num_matches = match_patterns(pattern, species_pattern)
        match = num_matches > 0
        return match


def generate_species_from_molecule_type(molecule_type: MoleculeType):
    name = molecule_type.name

    molecule_type_components = [component for name, component in molecule_type.components.items(
    ) for _ in range(molecule_type.components_counts[name])]
    components_iters = [
        components_gen(component) for component in molecule_type_components]

    for aux2 in product(*components_iters):
        molecule = Molecule(name, list(aux2))
        pattern = Pattern([molecule])
        yield pattern
