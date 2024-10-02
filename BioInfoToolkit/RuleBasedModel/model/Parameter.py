import ast


class Parameter:
    name: str
    expression: str
    value: int|float
    comment: str

    def __init__(self, name: str, expression: str) -> None:
        self.name = name
        self.expression = expression
        self.value = 0
        self.comment = ""

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
