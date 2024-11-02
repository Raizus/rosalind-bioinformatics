
from BioInfoToolkit.RuleBasedModel.utils.network_parsers import parse_observables_group

import numpy as np

class ObservablesGroup:
    name: str
    #(species_id => weight
    weighted_species: dict[int, int]

    def __init__(self, name: str, weighted_species: dict[int, int]) -> None:
        self.name = name
        self.weighted_species = weighted_species

    @classmethod
    def from_declaration(cls, declaration: str) -> "ObservablesGroup":
        parsed_group = parse_observables_group(declaration)
        label = parsed_group['label']
        weighted_species = parsed_group['weighted_species']
        group = ObservablesGroup(label, weighted_species)
        return group

    def compute_concentration(self, concentrations: dict[int, int] | dict[int, float]):
        concentration = sum(w*concentrations[sp_id]
                            for sp_id, w in self.weighted_species.items())
        return concentration

    def group_str(self) -> str:
        sorted_keys = sorted(self.weighted_species.keys())
        sp_w_gen = ((sp_id, self.weighted_species[sp_id]) for sp_id in sorted_keys)
        out = ','.join(f'{w}*{sp_id}' if w !=1 else f'{sp_id}'
                       for sp_id, w in sp_w_gen)
        return out

    def __repr__(self) -> str:
        out = f"{self.name} {self.group_str()}"
        return out

    def get_weights_array(self, n: int):
        weights = np.zeros(n, dtype=np.int64)

        for sp_id, w in self.weighted_species.items():
            weights[sp_id-1] = w

        return weights
