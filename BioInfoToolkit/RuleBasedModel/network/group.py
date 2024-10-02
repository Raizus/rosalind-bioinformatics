
from BioInfoToolkit.RuleBasedModel.network.parsers import parse_observables_group


class ObservablesGroup:
    name: str
    # (species_id, weight)
    weighted_species: list[tuple[int, int]]

    def __init__(self, name: str, weighted_species: list[tuple[int, int]]) -> None:
        self.name = name
        self.weighted_species = weighted_species

    @classmethod
    def from_declaration(cls, declaration: str) -> "ObservablesGroup":
        parsed_group = parse_observables_group(declaration)
        label = parsed_group['label']
        species = parsed_group['species']
        group = ObservablesGroup(label, species)
        return group

    def __repr__(self) -> str:
        out = f"{self.name} {','.join(f'{w}*{sp_id}' for sp_id, w in self.weighted_species)}"
        return out
