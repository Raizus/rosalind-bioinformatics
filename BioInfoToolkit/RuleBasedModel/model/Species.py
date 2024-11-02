
from itertools import product
from typing import Any, Generator
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.Compartment import Compartments
from BioInfoToolkit.RuleBasedModel.model.Component import components_gen
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_seed_species
from BioInfoToolkit.RuleBasedModel.model.Pattern import Molecule, Pattern, \
    match_pattern_specie, node_pattern_matching_func
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import SeedSpeciesDict


class Species:
    pattern: Pattern
    expression: str
    conc: float

    def __init__(self,
        pattern: Pattern,
        expression: str,
        conc: float = 0.0,
    ) -> None:
        self.pattern = pattern
        # we must change the pattern so that implicit molecule compartments are added if
        # the aggregate compartment is explicit
        if pattern.aggregate_compartment:
            for mol in pattern.molecules:
                if not mol.compartment:
                    mol.compartment = pattern.aggregate_compartment

        self.conc = conc
        self.expression = expression

    @classmethod
    def from_dict(
        cls,
        parsed: SeedSpeciesDict,
        molecules: dict[str, MoleculeType] | None
    ) -> "Species":
        pattern = Pattern.from_dict(parsed["pattern"], molecules)
        expression = parsed["expression"]
        specie = Species(pattern, expression)
        return specie

    @classmethod
    def from_declaration(cls,
                         declaration: str,
                         molecules: dict[str, MoleculeType] | None) -> "Species":
        parsed = parse_seed_species(declaration)
        specie = cls.from_dict(parsed, molecules)
        return specie

    def validate(self,
                 molecule_types: dict[str, MoleculeType],
                 compartments: Compartments) -> bool:
        """Checks if the specie is fully defined

        Args:
            molecule_types (dict[str, MoleculeType]): _description_

        Returns:
            bool: _description_
        """
        pattern = self.pattern
        if not pattern.is_connected():
            return False

        for mol in pattern.molecules:
            name = mol.name
            mol_type = molecule_types[name]

            # check if all components are fully defined
            if mol.components_counts != mol_type.components_counts:
                return False
            for comp in mol.components:
                if comp.states is None:
                    return False
                if not comp.is_stateless() and comp.state not in comp.states:
                    return False
                if comp.bond in ('?', '+'):
                    return False

        # when validating the species, the aggregate component could be implicit
        # in which case it should be set
        if all(mol.compartment for mol in pattern.molecules):
            comps = set(mol.compartment for mol in pattern.molecules if mol.compartment)
            aggregate_comp = compartments.get_aggregate(comps)
            if not pattern.aggregate_compartment:
                pattern.aggregate_compartment = aggregate_comp
            if aggregate_comp != pattern.aggregate_compartment:
                return False

        if not pattern.is_topologically_consistent(compartments):
            return False

        # for compartmentalized models, each molecule should have their compartment
        # explicitly declared, not just the aggregate
        if not compartments.is_compartmentalized():
            return True
        
        for mol in pattern.molecules:
            if not mol.compartment:
                return False
            if not compartments.has_compartment(mol.compartment):
                return False

        return True

    def __repr__(self) -> str:
        out = f"{str(self.pattern)} {self.expression}"
        return out

    def match_pattern(self, pattern: Pattern) -> bool:
        """Matched the species to a given pattern. To match, there must be a subgraph of the 
        species pattern graph that is isomorphic to the given pattern graph

        Args:
            pattern (Pattern): _description_

        Returns:
            bool: _description_
        """
        species_pattern = self.pattern
        num_matches = match_pattern_specie(pattern, species_pattern)
        match = num_matches > 0
        return match

    def iterate_matches(
            self,
            pattern: Pattern
        ) -> Generator[dict[tuple[int, int], tuple[int, int]], Any, None]:
        """Generates all subraph isomorphisms between a subgraph of the species graph and
        the graph of the given pattern

        Args:
            pattern (Pattern): _description_

        Yields:
            Generator[dict[tuple[int, int], tuple[int, int]], Any, None]: generator
            over the subgraph isomorphisms
        """
        matcher = nx.isomorphism.GraphMatcher(
            self.pattern.graph, pattern.graph, node_pattern_matching_func)

        yield from matcher.subgraph_isomorphisms_iter()


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
        cache (dict[str, int], optional): memoization cache. Maps str(pattern) to species_id. 
        Defaults to None.

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
