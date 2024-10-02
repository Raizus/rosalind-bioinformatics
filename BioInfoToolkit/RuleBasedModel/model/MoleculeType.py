from collections import Counter
from BioInfoToolkit.RuleBasedModel.model.Component import MoleculeTypeComponent, components_all_equal
from BioInfoToolkit.RuleBasedModel.model.Parsers import parse_molecule_type


class MoleculeType:
    name: str
    components_counts: Counter[str]
    components: dict[str, MoleculeTypeComponent]

    def __init__(self, name: str, components: list[MoleculeTypeComponent]) -> None:
        self.name = name
        self.components = {}

        # validate
        counts = Counter([comp.name for comp in components])
        unique_components: list[MoleculeTypeComponent] = []
        for comp_name in counts:
            comps = [comp for comp in components if comp.name == comp_name]
            if not components_all_equal(comps):
                raise ValueError("Components with the same name must have the same states.")
            unique_components.append(comps[0])

        self.components = {comp.name: comp for comp in unique_components}
        self.components_counts = counts

    @classmethod
    def from_declaration(cls, declaration: str) -> "MoleculeType":
        parsed = parse_molecule_type(declaration)
        if not parsed:
            raise ValueError(f"Invalid molecule declaration: {declaration}")

        name = parsed["name"]
        parsed_components = parsed['components']
        components = [MoleculeTypeComponent(p_comp['name'], p_comp['states'])
                      for p_comp in parsed_components]

        molecule_type = MoleculeType(name, components)
        return molecule_type

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
