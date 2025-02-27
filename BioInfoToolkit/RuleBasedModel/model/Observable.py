from BioInfoToolkit.RuleBasedModel.model.compartment import Compartments
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_observable
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern, match_pattern_specie
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import ObservableExpressionDict
from BioInfoToolkit.RuleBasedModel.utils.utls import apply_inequality


class ObservableElement:
    """Defines an Observable element, from which observables are formed.
    Pattern defines a pattern to compare with a given species.
    Sign and value can be used to find species of molecules that satisfy a given condition.
    For example, for an element defined as "R>20", the observable element will match with 
    all species that form a complex with more than 20 R molecules.

    Returns:
        ObservableElement: _description_
    """
    pattern: Pattern
    sign: str | None
    value: int | None

    def __init__(self,
                 pattern: Pattern,
                 sign: str | None,
                 value: int | None) -> None:
        self.pattern = pattern
        self.sign = sign
        self.value = value

    @classmethod
    def from_dict(cls, parsed: ObservableExpressionDict,
                  molecule_types: dict[str, MoleculeType] | None):
        parsed_pattern = parsed['pattern']
        sign = parsed['sign']
        value = parsed['value']
        pattern = Pattern.from_dict(parsed_pattern, molecule_types)

        element = ObservableElement(pattern, sign, value)
        return element

    def __repr__(self) -> str:
        out = f"{self.pattern}"
        if self.sign and self.value is not None:
            out += f"{self.sign}{self.value}"
        return out


class Observable:
    type: str
    label: str
    elements: list[ObservableElement]

    def __init__(self, obs_type: str, label: str,
                 elements: list[ObservableElement]) -> None:
        if obs_type not in ('Molecules', 'Species'):
            raise ValueError(f"type must be 'Molecules' or 'Species' (found {obs_type}).")
        self.type = obs_type
        self.label = label
        self.elements = elements

    def validate(
        self,
        molecule_types: dict[str, MoleculeType],
        compartments: Compartments
    ) -> bool:
        for element in self.elements:
            if not element.pattern.validate(molecule_types, compartments):
                return False

        return True

    @classmethod
    def from_declaration(cls, declaration: str, molecules: dict[str, MoleculeType]):
        parsed = parse_observable(declaration)
        if not parsed:
            raise ValueError(f"Invalid observable declaration: {declaration}")

        obs_type = parsed["type"]
        label = parsed["label"]
        parsed_expressions = parsed["expressions"]
        elements: list[ObservableElement] = []
        for parsed_expression in parsed_expressions:
            element = ObservableElement.from_dict(parsed_expression, molecules)
            elements.append(element)
        observable = Observable(obs_type, label, elements)
        return observable

    def match_species(self, species: Pattern) -> int:
        """Matches the observable object to a given species.
        A molecules-observable weights the concentrations of the species by 
        the number of matches to a pattern (or list of patterns).

        A species-observable is an unweighted sum of the concentrations of species
        matched to a pattern (or list of patterns).

        Args:
            species (Pattern): _description_

        Returns:
            int: If species-observable, returns 1 if there's a match, 0 otherwise.
            If molecules-observable, return 0 if there's no match, n otherwise, 
            where n is the number of all subgraph isomorphisms.
        """
        count_unique = not self.type == 'Species'
        total_matches = 0

        for element in self.elements:
            sign = element.sign
            value = element.value

            count_unique2 = count_unique or (bool(sign) and value is not None)
            num_matches = match_pattern_specie(element.pattern, species, count_unique2)

            if sign and value is not None:
                result = apply_inequality(num_matches, sign, value)
                if not result:
                    num_matches = 0
                elif not count_unique and result:
                    num_matches = 1

            total_matches += num_matches
        return total_matches

    def __repr__(self) -> str:
        out = f"{self.type}\t\t{self.label}\t{self.elements})"
        return out
