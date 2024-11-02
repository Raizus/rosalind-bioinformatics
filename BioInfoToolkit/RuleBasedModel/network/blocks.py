

from typing import OrderedDict
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.utils.utls import eval_expr, format_data_into_lines
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

    def evaluate_parameters(self):
        evaluated_params: OrderedDict[str, float | int] = OrderedDict()

        # evaluate parameters
        while len(evaluated_params) != len(self.items):
            new_eval = False

            for name, param in self.items.items():
                if name in evaluated_params:
                    continue
                value, _ = eval_expr(param.expression, evaluated_params)
                if not isinstance(value, int) and not isinstance(value, float):
                    raise TypeError(
                        f"value '{value}' must be of type int or float.")
                evaluated_params[name] = value
                new_eval = True

            if not new_eval:
                break

        # raise error if there are parameters that could not be evaluated
        if len(evaluated_params) != len(self.items):
            names = set(self.items.keys()) - set(evaluated_params.keys())
            names_str = ', '.join(names)
            msg = ("Could not evaluate all expressions. " +
                   f"Specifically could not evaluate parameters: {names_str}")
            raise ValueError(msg)

        self.evaluated_params = evaluated_params

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
            group_str = group.group_str()
            data.append((str(g_id), str(group.name), group_str))

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)
