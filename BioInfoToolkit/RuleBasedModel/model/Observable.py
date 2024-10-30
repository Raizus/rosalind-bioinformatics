from BioInfoToolkit.RuleBasedModel.model.Compartment import Compartments
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_observable
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern, match_pattern_specie
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import ObservableExpressionDict


class ObservableElement:
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
        count = not self.type == 'Species'
        total_matches = 0

        for element in self.elements:
            sign = element.sign
            value = element.value

            num_matches = match_pattern_specie(element.pattern, species, count,
                                               sign, value)

            total_matches += num_matches
        return total_matches

    def __repr__(self) -> str:
        out = f"{self.type}\t\t{self.label}\t{self.elements})"
        return out
