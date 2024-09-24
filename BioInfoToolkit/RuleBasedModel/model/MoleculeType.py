from collections import defaultdict
from BioInfoToolkit.RuleBasedModel.model.Component import MoleculeTypeComponent
from BioInfoToolkit.RuleBasedModel.model.Parsers import parse_molecule_type


class MoleculeType:
    name: str
    components_counts: defaultdict[str, int]
    components: dict[str, MoleculeTypeComponent]

    def __init__(self, declaration: str) -> None:
        parsed = parse_molecule_type(declaration)
        if not parsed:
            raise ValueError(f"Invalid molecule declaration: {declaration}")

        self.name = parsed["name"]
        self.components = dict()
        self.components_counts = defaultdict(int)

        parsed_components = parsed["components"]
        for parsed_component in parsed_components:
            comp_name = parsed_component["name"]
            states = parsed_component['states']
            component = MoleculeTypeComponent(comp_name, states)
            if comp_name in self.components:
                if self.components[comp_name] == component:
                    self.components_counts[comp_name] += 1
                else:
                    raise ValueError(
                        f"Component {component} has the same name as an already existing component but have different sets of allowed states {self.components[comp_name]}")
            else:
                self.components_counts[comp_name] += 1
                self.components[comp_name] = component

    def __repr__(self) -> str:
        out = f"{self.name}("
        for i, (comp_name, comp_count) in enumerate(self.components_counts.items()):
            component = self.components[comp_name]
            comp_str = str(component)
            comps_str = ','.join(comp_str for _ in range(comp_count))
            if i == 0:
                out += comps_str
            else:
                out += ',' + comps_str
        out += ')'
        return out
