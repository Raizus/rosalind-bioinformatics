
from itertools import product
from typing import Any, Generator
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.Component import components_gen
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_seed_species
from BioInfoToolkit.RuleBasedModel.model.Pattern import Molecule, Pattern, \
    match_pattern_specie, node_pattern_matching_func

class Species:
    compartment: str|None
    pattern: Pattern
    expression: str
    conc: float

    def __init__(self,
                 pattern: Pattern,
                 expression: str,
                 conc: float = 0.0,
                 compartment: str | None = None) -> None:
        # if not pattern.is_specie():
        #     raise ValueError(f"Reactant {pattern} must be a specie and not a pattern.")
        self.pattern = pattern
        self.conc = conc
        self.expression = expression
        self.compartment = compartment

    @classmethod
    def from_declaration(cls,
                         declaration: str,
                         molecules: dict[str, MoleculeType]):
        # parse string declaration
        parsed = parse_seed_species(declaration)
        if not parsed:
            raise ValueError(f"Invalid reactant declaration: {declaration}")

        pattern = Pattern.from_dict(parsed["pattern"], molecules)
        expression = parsed["expression"]
        specie = Species(pattern, expression)
        return specie

    def validate(self, molecule_types: dict[str, MoleculeType]) -> bool:
        """Checks if the specie is fully defined

        Args:
            molecule_types (dict[str, MoleculeType]): _description_

        Returns:
            bool: _description_
        """
        pattern = self.pattern
        for mol in pattern.molecules:
            if not pattern.is_connected():
                return False
            name = mol.name
            mol_type = molecule_types[name]
            if mol.components_counts != mol_type.components_counts:
                return False
            for comp in mol.components:
                if not comp.is_stateless() and comp.state not in comp.states:
                    return False
                if comp.bond in ('?', '+'):
                    return False

        return True

    def __repr__(self) -> str:
        out = ""
        if self.compartment:
            out = f"@{self.compartment}::"
        out += str(self.pattern)
        out += f" {self.expression}"
        return out

    def match_pattern(self, pattern: Pattern) -> bool:
        """Matched the species to a given pattern. To match, there must be a subgraph of the species pattern graph that is isomorphic to the given pattern graph

        Args:
            pattern (Pattern): _description_

        Returns:
            bool: _description_
        """
        species_pattern = self.pattern
        num_matches = match_pattern_specie(pattern, species_pattern)
        match = num_matches > 0
        return match

    def iterate_matches(self, pattern: Pattern
                        ) -> Generator[dict[tuple[int, int], tuple[int, int]], Any, None]:
        matcher = nx.isomorphism.GraphMatcher(
            self.pattern.graph, pattern.graph, node_pattern_matching_func)
        for mapping in matcher.subgraph_isomorphisms_iter():
            yield mapping

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


def find_species_match(pattern: Pattern, species: dict[int, Species]):
    """Returns the species_id that matches pattern exactly. Returns -1 if no match

    Args:
        pattern (Pattern): pattern to match
        species (dict[int, Species]): species dictionary

    Yields:
        int: sp_id (species_id)
    """
    for sp_id, specie in species.items():
        if specie.pattern == pattern:
            return sp_id
    return -1


def species_match_gen(pattern: Pattern,
                      species: dict[int, Species],
                      cache: dict[tuple[str,int], int] | None = None):
    """Generator of species_id that match the given pattern with memoization

    Args:
        pattern (Pattern): pattern to match
        species (dict[int, Species]): species dictionary
        cache (dict[str, int], optional): memoization cache. Defaults to None.

    Yields:
        int: sp_id (species_id)
    """
    if cache is None:
        cache = {}

    pattern_str = str(pattern)

    # Check if the result is already cached
    if pattern_str in cache:
        sp_id = cache[pattern_str]
        if sp_id is not None:
            yield sp_id
        return

    for sp_id, specie in species.items():
        key = (pattern_str, sp_id)
        if key in cache:
            yield cache[key]
        else:
            match = specie.match_pattern(pattern)
            if match:
                cache[(pattern_str, sp_id)] = sp_id
                yield sp_id
