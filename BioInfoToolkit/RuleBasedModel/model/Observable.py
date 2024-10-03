from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_observable
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern, match_pattern_specie


class Observable:
    type: str
    label: str
    patterns: list[Pattern]

    def __init__(self, obs_type: str, label: str, patterns: list[Pattern]) -> None:
        if obs_type not in ('Molecules', 'Species'):
            raise ValueError(f"type must be 'Molecules' or 'Species' (found {obs_type}).")
        self.type = obs_type
        self.label = label
        self.patterns = patterns

    def validate(self, molecule_types: dict[str, MoleculeType]) -> bool:
        for pattern in self.patterns:
            if not pattern.validate(molecule_types):
                return False

        return True

    @classmethod
    def from_declaration(cls, declaration: str, molecules: dict[str, MoleculeType]):
        parsed = parse_observable(declaration)
        if not parsed:
            raise ValueError(f"Invalid observable declaration: {declaration}")

        obs_type = parsed["type"]
        label = parsed["label"]
        parsed_patterns = parsed["patterns"]
        patterns: list[Pattern] = []
        for parsed_pattern in parsed_patterns:
            pattern = Pattern.from_dict(parsed_pattern, molecules)
            patterns.append(pattern)
        observable = Observable(obs_type, label, patterns)
        return observable

    def match_species(self, species: Pattern) -> int:
        count = not self.type == 'Species'
        num_matches = 0
        for pattern in self.patterns:
            num_matches += match_pattern_specie(pattern, species, count)
        return num_matches

    def __repr__(self) -> str:
        out = f"{self.type}\t\t{self.label}\t{self.patterns})"
        return out
