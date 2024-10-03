from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_compartment
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import CompartmentDict


class Compartment:
    name: str
    dimensions: int
    volume: str
    enclosing_compartment: str | None
    comment: str

    def __init__(self, name: str, dimensions: int,
                 volume: str, enclosing_compartment: str | None,
                 comment: str) -> None:
        self.name = name
        self.dimensions = dimensions
        self.volume = volume
        self.enclosing_compartment = enclosing_compartment
        self.comment = comment

    @classmethod
    def from_dict(cls, parsed: CompartmentDict) -> "Compartment":
        name = parsed["name"]
        dimensions = parsed["dimensions"]
        volume = parsed["volume"]
        enclosing_compartment = parsed["enclosing_compartment"]
        comment = parsed["comment"]

        compartment = Compartment(name, dimensions, volume, enclosing_compartment, comment)
        return compartment

    @classmethod
    def from_declaration(cls, declaration: str) -> "Compartment":
        parsed = parse_compartment(declaration)
        compartment = cls.from_dict(parsed)
        return compartment

    def is_enclosed(self) -> bool:
        return bool(self.enclosing_compartment)

    def __repr__(self) -> str:
        out = f"{self.name} {self.dimensions} {self.volume} {self.enclosing_compartment}"
        return out
