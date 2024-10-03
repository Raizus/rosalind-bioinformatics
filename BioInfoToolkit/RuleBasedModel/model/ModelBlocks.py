

from typing import OrderedDict

from BioInfoToolkit.RuleBasedModel.model.Compartment import Compartment
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
        """Returns the string of this model block.

        Returns:
            str: output string
        """
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
        """Returns the string of this model block.

        Returns:
            str: output string
        """
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str,str,str]] = []
        for _, observable in self.items.items():
            data.append((observable.type, observable.label, str(observable.patterns)))

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
        """Returns the string of this model block.

        Returns:
            str: output string
        """
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
        """Returns the string of this model block.

        Returns:
            str: output string
        """
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
        """Returns the string of this model block.

        Returns:
            str: output string
        """
        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str]] = []
        for _, reaction in self.items.items():
            name = f"{reaction.name}:" if reaction.name else ''
            reactants = reaction.reactants
            products = reaction.products
            reactants_str = ' + '.join(str(reagent) for reagent in reactants)
            products_str = ' + '.join(str(reagent) for reagent in products)
            arrow = '<->' if reaction.is_bidirectional() else '->'
            reaction_str = f"{reactants_str} {arrow} {products_str}"

            rate_f = reaction.forward_rate
            rate_r = reaction.reverse_rate
            rates_str = f"{rate_f}"
            if rate_r:
                rates_str = rates_str + f", {rate_r}"
            data_line = (name, reaction_str, rates_str)
            data.append(data_line)

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class CompartmentsBlock(ModelBlock):
    name = 'compartments'
    items: OrderedDict[str, Compartment]

    def __init__(self) -> None:
        super().__init__()
        self.items = OrderedDict()

    def add_compartment(self, compartment: Compartment):
        name = compartment.name
        dimensions = compartment.dimensions
        enclosing_compartment_name = compartment.enclosing_compartment

        if name in self.items:
            raise ValueError(f"Compartment with {name} already declared.")

        if len(self.items) == 0:
            if compartment.is_enclosed():
                msg = ( "The first compartment must be the external one, " +
                        "enclosing_compartment must be None.")
                raise TypeError(msg)
        else:
            if not compartment.is_enclosed() or (enclosing_compartment_name not in self.items):
                msg = (f"The compartment '{name}' must be enclosed by another" +
                       "that has been declared. " +
                       f"Compartment '{enclosing_compartment_name}' not declared.")
                raise ValueError(msg)

            # this compartment is a surface and must be enclosed by a volumne
            enclosing_compartment = self.items[enclosing_compartment_name]
            if dimensions == 2:
                if enclosing_compartment.dimensions != 3:
                    msg = f"The compartment '{name}' is a surface and must be enclosed by a volume."
                    raise ValueError(msg)

            else: # this compartment is a volume
                if enclosing_compartment.dimensions != 2:
                    msg = f"The compartment '{name}' is a volume and must be enclosed by a volume."
                    raise ValueError(msg)

        # TODO build compartment tree, and check that every surface encloses exactly one volume

        self.items[name] = compartment

    def gen_string(self) -> str:
        """Returns the string of this model block.

        Returns:
            str: output string
        """
        if len(self.items) == 0:
            return ''

        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str, str, str]] = []
        for _, compartment in self.items.items():
            name = compartment.name
            dimensions = compartment.dimensions
            volume = compartment.volume
            enclosing_compartment = compartment.enclosing_compartment
            enclosing_compartment = enclosing_compartment if enclosing_compartment else ''
            comment = f"# {compartment.comment}" if compartment.comment else ''
            data_line = (name, str(dimensions), volume, enclosing_compartment, comment)
            data.append(data_line)

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)
