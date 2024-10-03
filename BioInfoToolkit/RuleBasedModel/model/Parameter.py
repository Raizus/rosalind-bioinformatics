import ast

from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_parameter
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import ParameterDict


class Parameter:
    name: str
    expression: str
    value: int|float
    comment: str

    def __init__(self, name: str, expression: str, comment: str = '') -> None:
        self.name = name
        self.expression = expression
        self.value = 0
        self.comment = comment

    @classmethod
    def from_dict(cls, parsed: ParameterDict) -> "Parameter":
        name = parsed['name']
        expression = parsed['expression']
        comment = parsed['comment']
        parameter = Parameter(name, expression, comment)
        return parameter

    @classmethod
    def from_declaration(cls, declaration: str) -> "Parameter":
        parsed = parse_parameter(declaration)
        parameter = cls.from_dict(parsed)
        return parameter

    def __repr__(self) -> str:
        out = f"{self.name} {self.expression}"
        return out

    def eval(self) -> bool:
        expr = self.expression
        try:
            value = ast.literal_eval(expr)
            self.value = value
            self.comment = "Constant"
            return True
        except Exception:
            return False
            # raise ValueError("Some error.") from exc
            # value = eval(declaration, {"__builtins__": None}, dict_aux)

    def copy(self) -> "Parameter":
        new_param = Parameter(self.name, self.expression)
        return new_param
