from typing import Any, TypedDict
import pyparsing as pp

# Define the grammar for parsing
pp.ParserElement.enablePackrat()

# Define a word starting with a letter
NAME_EXPRESSION = pp.Word(pp.alphas, pp.alphanums + "_")
NUMS = pp.Word(pp.nums)

# Define state and bond (alphanumeric, can start with number)
state = NAME_EXPRESSION ^ NUMS
BOND_NAME = NUMS ^ pp.Literal('+') ^ pp.Literal('?')
BOND_EXPRESSION = pp.Optional(pp.Literal(
    '!') + BOND_NAME('bond')).leaveWhitespace()

# Define component format
COMPONENT_PARSER = pp.Group(
    NAME_EXPRESSION('name') +
    pp.Optional(pp.Literal('~') + state('state')).leaveWhitespace() +
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

variable = pp.Word(pp.alphas, pp.alphanums+'_')

# Define numeric literals
number = pp.Word("+-" + pp.nums + ".", pp.nums + ".")

# Operators: +, -, *, /, and ^
plus = pp.oneOf("+ -")
mult = pp.oneOf("* /")
exp = pp.Literal("^")

# Forward declaration for handling nested parentheses
EXPRESSION_PATTERN = pp.Forward()

# Parentheses
lparen = pp.Suppress("(")
rparen = pp.Suppress(")")

# Define an expression that can be a number, variable, or parenthesized expression
operand = number | variable | (lparen + EXPRESSION_PATTERN + rparen)

# Define the expression grammar using infix notation
EXPRESSION_PATTERN <<= pp.infixNotation(
    operand,
    [
        (exp, 2, pp.opAssoc.RIGHT),  # Exponentiation
        (mult, 2, pp.opAssoc.LEFT),  # Multiplication and division
        (plus, 2, pp.opAssoc.LEFT),  # Addition and subtraction
    ]
)


class MoleculeTypeComponentDict(TypedDict):
    name: str
    states: set[str]


class MoleculeTypeDict(TypedDict):
    name: str
    components: list[MoleculeTypeComponentDict]


class ComponentDict(TypedDict):
    name: str
    state: str
    bond: str


class MoleculeDict(TypedDict):
    name: str
    components: list[ComponentDict]


class ReactionDict(TypedDict):
    name: str
    reactants: list[list[MoleculeDict]]
    products: list[list[MoleculeDict]]
    forward_rate: str
    reverse_rate: str | None


class ObservableDict(TypedDict):
    type: str
    label: str
    pattern: list[MoleculeDict]


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


def parse_molecule(declaration: str):
    molecule_parser = MOLECULE_PARSER.leaveWhitespace() + pp.StringEnd()
    try:
        parsed = molecule_parser.parseString(declaration)
    except pp.ParseException as ee:
        raise ValueError(
            f"Reactant {declaration} not declared correctly") from ee

    result = parsed_molecule_to_dict(parsed)
    return result


def parse_pattern(declaration: str) -> list[MoleculeDict]:
    try:
        parsed = PATTERN_PARSER.parseString(declaration)
    except pp.ParseException as ee:
        raise ValueError(
            f"Complex {declaration} not declared correctly.") from ee

    if not isinstance(parsed, pp.ParseResults):
        raise TypeError(f"{parsed} must be of type ParseResults.")

    molecules = parsed_pattern_to_dict_list(parsed)
    return molecules


def parsed_reactants_to_list(parsed: pp.ParseResults | Any):
    if not isinstance(parsed, pp.ParseResults):
        raise TypeError(f"{parsed} must be of type ParseResults.")

    patterns: list[list[MoleculeDict]] = []
    for parsed_pattern in parsed:
        if not isinstance(parsed_pattern, pp.ParseResults):
            raise TypeError(f"{parsed_pattern} must be of type ParseResults.")

        parts: list[MoleculeDict] = []
        for parsed_molecule in parsed_pattern:
            reactant = parsed_molecule_to_dict(parsed_molecule)
            parts.append(reactant)

        patterns.append(parts)

    return patterns


def parse_reactants_sum(string: str):
    reagent_pattern = pp.Group(
        COMPLEX_PARSER | pp.Group(MOLECULE_PARSER))
    reagents_pattern = pp.delimitedList(reagent_pattern, '+')

    left_side_reactants = reagents_pattern('reactants')
    parsed = left_side_reactants.parseString(string)

    reactants: list[list[MoleculeDict]
                    ] = parsed_reactants_to_list(parsed.reactants)

    return reactants



def parsed_reaction_to_reaction_dict(parsed: pp.ParseResults):
    name = parsed.name
    forward_rate = parsed.forward_rate
    reverse_rate = parsed.get('reverse_rate', None)

    if not isinstance(name, str):
        raise ValueError("Name must be a string.")

    if not isinstance(forward_rate, str):
        raise ValueError("forward_rate must be a string.")

    if reverse_rate is not None and not isinstance(reverse_rate, str):
        raise ValueError("reverse_rate must be a string or none.")

    left_side = parsed.left_side_reactants
    right_side = parsed.right_side_reactants

    left_reactants = parsed_reactants_to_list(left_side)
    right_reactants = parsed_reactants_to_list(right_side)

    result: ReactionDict = {
        "name": name,
        "reactants": left_reactants,
        "products": right_reactants,
        "forward_rate": forward_rate,
        "reverse_rate": reverse_rate
    }

    return result


def parse_reaction(reaction_str: str):
    """Parses a biochemical reaction and verifies the correctness of complexes and bonds."""

    reagent_pattern = pp.Group(
        COMPLEX_PARSER | pp.Group(MOLECULE_PARSER))
    reagents_pattern = pp.delimitedList(reagent_pattern, '+')

    left_side_reactants = reagents_pattern('left_side_reactants')
    right_side_reactants = reagents_pattern('right_side_reactants')

    arrow_unidirectional = pp.Suppress('->')
    arrow_bidirectional = pp.Suppress('<->')

    expression_parser = pp.Combine(EXPRESSION_PATTERN)

    forward_rate_pattern = expression_parser('forward_rate')
    reverse_rate_pattern = expression_parser('reverse_rate')

    unidirectional_pattern = (NAME_EXPRESSION('name') + pp.Suppress(':')).leaveWhitespace() + \
        left_side_reactants + \
        arrow_unidirectional + \
        right_side_reactants + \
        forward_rate_pattern

    bidirectional_pattern = (NAME_EXPRESSION('name') + pp.Literal(':')).leaveWhitespace() + \
        left_side_reactants + \
        arrow_bidirectional + \
        right_side_reactants + \
        forward_rate_pattern + \
        pp.Literal(',') + \
        reverse_rate_pattern

    rule_parser = unidirectional_pattern | bidirectional_pattern

    try:
        parsed = rule_parser.parseString(reaction_str)
        reaction = parsed_reaction_to_reaction_dict(parsed)
        return reaction
    except (pp.ParseException, pp.ParseBaseException) as ee:
        raise ValueError(
            f"Reaction {reaction_str} not declared correctly.") from ee


def parse_molecule_type(declaration: str):
    """A valid declaration has the format: name() or name(parameters).
    Parameters are comma separated, without spaces after the comma, and can be a word specifying a binding site, or a word followed by one or more tilde and a string specifying a state and possible state values. For example:
        - T(l,Phos~U~P) means T has a binding site with L, and the can be either phosphorylated or unphosphorylated.
        - or T(l,r,Phos~U~P,Meth~A~B~C) same has above, plus a binding side with R and a methilation state with can be either A, B or C
    Returns a dictionary with name, binding_sites, and a dictionary of state dict[str, set[str]], if the name is valid    

    Args:
        declaration (str): _description_

    Returns:
        _type_: _description_
    """

    # Define the format for components enclosed in parentheses
    molecule_type_state = pp.Suppress('~') + state
    molecule_type_states = pp.Group(pp.Optional(molecule_type_state[2, ...]))
    molecule_type_component = \
        (NAME_EXPRESSION('name') + molecule_type_states('states')).leaveWhitespace()

    molecule_type_components = pp.delimitedList(
        pp.Group(molecule_type_component)).leaveWhitespace()

    molecule_type_parser = (
        NAME_EXPRESSION('molecule_name') +
        pp.Literal('(').leaveWhitespace() +
        pp.Optional(molecule_type_components)('components') +
        pp.Literal(')').leaveWhitespace()
    )

    try:
        parsed = molecule_type_parser.parseString(declaration)
    except pp.ParseException as e:
        raise ValueError(
            f"Molecule {declaration} not declared correctly.") from e

    name = parsed.molecule_name
    if not isinstance(name, str):
        raise ValueError("Molecule name must be a string.")

    components: list[MoleculeTypeComponentDict] = []
    parsed_components = parsed.components
    if isinstance(parsed_components, pp.ParseResults):
        for parsed_comp in parsed_components:
            states = set(parsed_comp.states.as_list())
            comp: MoleculeTypeComponentDict = {
                'name': parsed_comp.name,
                'states': states,
            }
            components.append(comp)

    result: MoleculeTypeDict = {
        'name': name,
        'components': components
    }
    return result


def parse_observable(declaration: str) -> ObservableDict:
    type_pattern = pp.Literal('Molecules') | pp.Literal('Species')
    label_pattern = pp.Word(pp.alphas, pp.alphanums + '_')

    reactant_pattern = pp.Group(
        COMPLEX_PARSER | pp.Group(MOLECULE_PARSER))

    observable_parser = type_pattern(
        'type') + label_pattern('label') + reactant_pattern('pattern')
    # + pp.StringEnd()

    parsed = observable_parser.parseString(declaration)
    obs_type = parsed.type
    label = parsed.label

    if not isinstance(obs_type, str):
        raise TypeError(f"type {obs_type} must be a string")

    if not isinstance(label, str):
        raise TypeError(f"label {label} must be a string")

    parsed_reactant = parsed.pattern
    pattern = parsed_pattern_to_dict_list(parsed_reactant)

    result: ObservableDict = {
        'type': obs_type,
        'label': label,
        'pattern': pattern,
    }

    return result
