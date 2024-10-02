

from typing import OrderedDict

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Observable import Observable
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule
from BioInfoToolkit.RuleBasedModel.model.Species import Species
from BioInfoToolkit.RuleBasedModel.utils.utls import format_data_into_lines


class ModelBlock:
    name: str
    items: OrderedDict


class MoleculeTypesBlock(ModelBlock):
    name = "molecule types"
    items: OrderedDict[str, MoleculeType]

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()

    def add_molecule_type(self, mol_type: MoleculeType):
        if mol_type.name in self.items:
            raise ValueError(
                f"Molecule with name '{mol_type.name}' already declared.")

        self.items[mol_type.name] = mol_type

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        for _, mol_type in self.items.items():
            lines.append(f"\t{mol_type}")

        lines.append(f"end {self.name}\n")
        return '\n'.join(lines)


class ObservablesBlock(ModelBlock):
    name = "observables"
    items: OrderedDict[str, Observable]

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()

    def add_observable(self, observable: Observable):
        if observable.label in self.items:
            raise ValueError(
                f"Observable with label '{observable.label}' already declared.")
        self.items[observable.label] = observable

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str,str,str]] = []
        for _, observable in self.items.items():
            data.append((observable.type, observable.label, str(observable.pattern)))

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class ParametersBlock(ModelBlock):
    name = 'parameters'
    items: OrderedDict[str, Parameter]

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()

    def add_parameter(self, parameter: Parameter):
        self.items[parameter.name] = parameter
        # value = eval(declaration, {"__builtins__": None}, dict_aux)

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str]] = []
        for name, param in self.items.items():
            data.append((name, param.expression))

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class SeedSpeciesBlock(ModelBlock):
    name = 'species'
    items: OrderedDict[int, Species]
    id_count: int

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()
        self.id_count = 1

    def add_species(self, species: Species):
        self.items[self.id_count] = species
        self.id_count += 1

    def validate_species(self, molecule_types: dict[str, MoleculeType]) -> bool:
        for _, specie in self.items.items():
            valid = specie.validate(molecule_types)
            if not valid:
                return False
        return True

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str]] = []
        for _, seed_species in self.items.items():
            data.append((str(seed_species.pattern), seed_species.expression))

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class ReactionRulesBlock(ModelBlock):
    name = "reaction rules"
    items: OrderedDict[int, ReactionRule]
    count_id: int

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()
        self.count_id = 1

    def add_rule(self, rule: ReactionRule):
        # if rule in self.items:
        #     raise ValueError(
        #         f"Observable with label '{observable.label}' already declared.")
        self.items[self.count_id] = rule
        self.count_id += 1

    def gen_string(self):
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        for _, rule in self.items.items():
            lines.append(f"\t{rule}")

        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)
