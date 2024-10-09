
from BioInfoToolkit.RuleBasedModel.utils.network_parsers import parse_observables_group

import numpy as np

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

    def compute_concentration(self, concentrations: dict[int, int] | dict[int, float]):
        concentration = sum(w*concentrations[sp_id]
                            for sp_id, w in self.weighted_species)
        return concentration

    def __repr__(self) -> str:
        out = f"{self.name} {','.join(f'{w}*{sp_id}' for sp_id, w in self.weighted_species)}"
        return out

    def get_weights_array(self, n: int):
        weights = np.zeros(n, dtype=np.int64)

        for sp_id, w in self.weighted_species:
            weights[sp_id-1] = w

        return weights
