from typing import Any, TypedDict

import pyparsing as pp


class MoleculeTypeComponentDict(TypedDict):
    """Typed dict for molecule type component"""
    name: str
    states: set[str]


class MoleculeTypeDict(TypedDict):
    """Typed dict for molecule type"""
    name: str
    components: list[MoleculeTypeComponentDict]


class ComponentDict(TypedDict):
    """Typed dict for component"""
    name: str
    state: str
    bond: str


class MoleculeDict(TypedDict):
    """Typed dict for molecule"""
    name: str
    components: list[ComponentDict]


class ReactionRuleDict(TypedDict):
    """Typed dict for reaction rule"""
    name: str
    reactants: list[list[MoleculeDict]]
    products: list[list[MoleculeDict]]
    forward_rate: str
    reverse_rate: str | None


class ReactionDict(TypedDict):
    """Typed dict for reaction"""
    id: int
    reactants: list[int]
    products: list[int]
    rate: str
    comment: str


class ObservableDict(TypedDict):
    """Typed dict for observable"""
    type: str
    label: str
    pattern: list[MoleculeDict]


class ObservablesGroupDict(TypedDict):
    """Typed dict for observable group"""
    id: int
    label: str
    species: list[tuple[int,int]]


class ParameterDict(TypedDict):
    """Typed dict for parameter"""
    name: str
    expression: str
    comment: str


class NetworkParameterDict(ParameterDict):
    """Typed dict for network parameter"""
    id: int


class SeedSpeciesDict(TypedDict):
    """Typed dict for seed species"""
    pattern: list[MoleculeDict]
    expression: str


class NetworkSeedSpeciesDict(SeedSpeciesDict):
    """Typed dict for seed species"""
    id: int


# Define the grammar for parsing
pp.ParserElement.enablePackrat()

# Define a word starting with a letter
NAME_EXPRESSION = pp.Word(pp.alphas, pp.alphanums + "_")
NUMS = pp.Word(pp.nums)
# Define state and bond (alphanumeric, can start with number)
STATE = NAME_EXPRESSION ^ NUMS
BOND_NAME = NUMS ^ pp.Literal('+') ^ pp.Literal('?')
BOND_EXPRESSION = pp.Optional(pp.Literal(
    '!') + BOND_NAME('bond')).leaveWhitespace()
# Define component format
COMPONENT_PARSER = pp.Group(
    NAME_EXPRESSION('name') +
    pp.Optional(pp.Literal('~') + STATE('state')).leaveWhitespace() +
    pp.Optional(pp.Literal('!') + BOND_NAME('bond')).leaveWhitespace()
)
# Define the format for components enclosed in parentheses
COMPONENTS_PARSER = pp.delimitedList(COMPONENT_PARSER, ',').leaveWhitespace()
COMBINED_COMPONENTS_PARSER = pp.delimitedList(
    COMPONENT_PARSER, ',', True).leaveWhitespace()
# Define the overall expression format
MOLECULE_PARSER = (NAME_EXPRESSION('molecule_name') +
                   pp.Literal('(').leaveWhitespace() +
                   pp.Optional(COMPONENTS_PARSER)('components') +
                   pp.Literal(')').leaveWhitespace())
COMBINED_MOLECULE_PARSER = pp.Combine(
    NAME_EXPRESSION('molecule_name') +
    pp.Literal('(').leaveWhitespace() +
    pp.Optional(COMBINED_COMPONENTS_PARSER)('components') +
    pp.Literal(')').leaveWhitespace()
)
COMPLEX_PARSER = pp.delimitedList(pp.Group(MOLECULE_PARSER), '.', min=2)
COMBINED_COMPLEX_PARSER = pp.delimitedList(
    pp.Group(COMBINED_MOLECULE_PARSER), '.', True, min=2)
# COMPLEX_PARSER = pp.Group(MOLECULE_PARSER) + \
#     pp.OneOrMore(pp.Suppress('.').leaveWhitespace() + \
#     pp.Group(MOLECULE_PARSER).leaveWhitespace())

PATTERN_PARSER = COMPLEX_PARSER | pp.Group(MOLECULE_PARSER)
LABEL_PARSER = pp.Word(pp.alphas, pp.alphanums+'_')

# Define numeric literals
sign = pp.Optional(pp.oneOf("+ -"))
integer_part = pp.Word(pp.nums)
decimal_part = pp.Optional(pp.Literal('.') + pp.Word(pp.nums))
# Scientific notation component
exp_part = pp.Optional(pp.CaselessLiteral(
    'e') + pp.Optional(pp.oneOf("+ -")) + pp.Word(pp.nums))

# Define the complete number parser
NUMBER_PARSER = pp.Combine(sign + integer_part + decimal_part + exp_part)

# Operators: +, -, *, /, and ^
PLUS_PARSER = pp.oneOf("+ -")
MULT_PARSER = pp.oneOf("* /")
EXP_PARSER = pp.Literal("^")

# Forward declaration for handling nested parentheses
EXPRESSION_PARSER = pp.Forward()

# Parentheses
LPAREN_PARSER = pp.Suppress("(")
RPAREN_PARSER = pp.Suppress(")")

# Define an expression that can be a number, variable, or parenthesized expression
OPERAND_PARSER = NUMBER_PARSER | LABEL_PARSER | (
    LPAREN_PARSER + EXPRESSION_PARSER + RPAREN_PARSER)

# Define the expression grammar using infix notation
EXPRESSION_PARSER <<= pp.infixNotation(
    OPERAND_PARSER,
    [
        (EXP_PARSER, 2, pp.opAssoc.RIGHT),  # Exponentiation
        (MULT_PARSER, 2, pp.opAssoc.LEFT),  # Multiplication and division
        (PLUS_PARSER, 2, pp.opAssoc.LEFT),  # Addition and subtraction
    ]
)

TEXT_PARSER = pp.Combine(pp.ZeroOrMore(pp.Word(pp.printables)))
COMMENT_PARSER = (pp.Suppress("#") + TEXT_PARSER('comment') +
                  (pp.LineEnd() | pp.StringEnd()))


def parsed_molecule_to_dict(parsed: pp.ParseResults):
    name = parsed.molecule_name
    if not isinstance(name, str):
        raise ValueError("Molecule name must be a string.")

    components: list[ComponentDict] = []
    parsed_components = parsed.components
    if isinstance(parsed_components, pp.ParseResults):
        for parsed_comp in parsed_components:
            comp: ComponentDict = {
                'name': parsed_comp.name,
                'state': parsed_comp.get('state', ''),
                'bond': parsed_comp.get('bond', '')
            }
            components.append(comp)

    result: MoleculeDict = {
        'name': name,
        'components': components
    }
    return result


def parsed_pattern_to_dict_list(parsed: pp.ParseResults | Any):
    if not isinstance(parsed, pp.ParseResults):
        raise TypeError(f"{parsed} must be of type ParseResults.")

    parts: list[MoleculeDict] = []
    for parsed_molecule in parsed:
        reactant = parsed_molecule_to_dict(parsed_molecule)
        parts.append(reactant)
    return parts


def parsed_parameter_to_parameter_dict(parsed: pp.ParseResults) -> ParameterDict:
    name = parsed.name
    if not isinstance(name, str):
        raise ValueError("Parameter name must be a string.")

    expression = parsed.expression
    if not isinstance(expression, str):
        raise ValueError("Parameter expression must be a string.")

    comment = parsed.get('comment', '')
    if not isinstance(comment, str):
        raise TypeError(f"type {comment} must be a string")

    result: ParameterDict = {
        'name': name,
        'expression': expression,
        'comment': comment
    }
    return result


def parsed_seed_species_to_seed_species_dict(parsed: pp.ParseResults) -> SeedSpeciesDict:
    parsed_pattern = parsed.pattern
    pattern = parsed_pattern_to_dict_list(parsed_pattern)

    expression = parsed.expression
    if not isinstance(expression, str):
        raise ValueError("See species expression must be a string.")

    result: SeedSpeciesDict = {
        'pattern': pattern,
        'expression': expression,
    }
    return result
