from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Parsers import parse_observable
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern, match_patterns


class Observable:
    type: str
    label: str
    pattern: Pattern

    def __init__(self, obs_type: str, label: str, pattern: Pattern) -> None:
        if obs_type != 'Molecules' and obs_type != 'Species':
            raise ValueError(f"type must be 'Molecules' or 'Species' (found {obs_type}).")
        self.type = obs_type
        self.label = label
        self.pattern = pattern

    @classmethod
    def from_declaration(cls, declaration: str, molecules: dict[str, MoleculeType]):
        parsed = parse_observable(declaration)
        if not parsed:
            raise ValueError(f"Invalid observable declaration: {declaration}")

        obs_type = parsed["type"]
        label = parsed["label"]
        parsed_pattern = parsed["pattern"]
        pattern = Pattern.from_dict(parsed_pattern, molecules)
        observable = Observable(obs_type, label, pattern)
        return observable

    def match_species(self, species: Pattern):
        count = not self.type == 'Species'
        num_matches = match_patterns(self.pattern, species, count)
        return num_matches
