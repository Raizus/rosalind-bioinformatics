

from typing import OrderedDict
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.utils.utls import format_data_into_lines
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup


class NetworkBlock:
    name: str

class ParametersBlock(NetworkBlock):
    name = "parameters"
    items: OrderedDict[str, Parameter]
    evaluated_params: OrderedDict[str, float | int]

    def __init__(self) -> None:
        self.items = OrderedDict()
        self.evaluated_params = OrderedDict()

    def add_parameter(self, parameter: Parameter):
        self.items[parameter.name] = parameter

    def gen_string(self):
        if len(self.items) == 0:
            return ''

        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str, str]] = []
        for i, (name, param) in enumerate(self.items.items()):
            comment = f"# {param.comment}" if param.comment else ''
            data_line = (str(i), name, param.expression, comment)
            data.append(data_line)

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)


class GroupsBlock(NetworkBlock):
    name = "groups"
    items: OrderedDict[int, ObservablesGroup]
    count_id: int

    def __init__(self) -> None:
        self.items = OrderedDict()
        self.count_id = 0

    def add_group(self, group: ObservablesGroup):
        self.items[self.count_id] = group
        self.count_id += 1

    def gen_string(self):
        if len(self.items) == 0:
            return ''

        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str]] = []
        for g_id, group in self.items.items():
            group_str = ','.join(f'{w}*{sp_id}' if w != 1 else f'{sp_id}' for sp_id,
                                 w in group.weighted_species)
            data.append((str(g_id), str(group.name), group_str))

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)
