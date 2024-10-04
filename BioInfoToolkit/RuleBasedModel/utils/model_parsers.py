from typing import Any
import pyparsing as pp

from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import COMMENT_PARSER, NUMS, UNSIGNED_NUMBER_PARSER, VARIABLE_PARSER, CompartmentDict, GenerateNetworkDict, MoleculeTypeComponentDict, \
    MoleculeTypeDict, MoleculeDict, ObservableExpressionDict, ParsingError, ReactionRuleDict, ObservableDict, \
    ParameterDict, SeedSpeciesDict, NAME_EXPRESSION, STATE, MOLECULE_PARSER, \
    COMPLEX_PARSER, PATTERN_PARSER, EXPRESSION_PARSER, parsed_compartment_to_compartment_dict, parsed_parameter_to_parameter_dict,parsed_seed_species_to_seed_species_dict
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import parsed_pattern_to_dict_list
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import parsed_molecule_to_dict

# Define the grammar for parsing
pp.ParserElement.enablePackrat()

def parsed_observable_to_dict(parsed: pp.ParseResults) -> ObservableDict:
    obs_type = parsed.type
    label = parsed.label

    if not isinstance(obs_type, str):
        raise TypeError(f"type {obs_type} must be a string")

    if not isinstance(label, str):
        raise TypeError(f"label {label} must be a string")

    parsed_expressions = parsed.expressions
    if not isinstance(parsed_expressions, pp.ParseResults):
        raise TypeError(f"{parsed_expressions} must be of type ParseResults.")

    expressions: list[ObservableExpressionDict] = []
    for parsed_exp in parsed_expressions:
        parsed_pattern = parsed_exp.pattern
        sign: str | None = parsed_exp.get('sign', None)
        parsed_value = parsed_exp.get('value', None)

        value: int | None = None
        if isinstance(parsed_value, str):
            value = int(parsed_value)

        pattern = parsed_pattern_to_dict_list(parsed_pattern)        
        obs_expr: ObservableExpressionDict = {
            'pattern': pattern,
            'sign': sign,
            'value': value
        }

        expressions.append(obs_expr)


    result: ObservableDict = {
        'type': obs_type,
        'label': label,
        'expressions': expressions
    }

    return result

def parse_molecule(declaration: str):
    molecule_parser = MOLECULE_PARSER.leaveWhitespace() + pp.StringEnd()
    try:
        parsed = molecule_parser.parseString(declaration)
    except pp.ParseException as ee:
        raise ParsingError(
            f"Reactant {declaration} not declared correctly") from ee

    result = parsed_molecule_to_dict(parsed)
    return result


def parse_pattern(declaration: str) -> list[MoleculeDict]:
    try:
        parsed = PATTERN_PARSER.parseString(declaration)
    except pp.ParseException as ee:
        raise ParsingError(f"Complex {declaration} not declared correctly.") from ee

    if not isinstance(parsed, pp.ParseResults):
        raise TypeError(f"{parsed} must be of type ParseResults.")

    molecules = parsed_pattern_to_dict_list(parsed)
    return molecules


def parsed_reactants_to_list(parsed: pp.ParseResults | Any) -> list[list[MoleculeDict]]:
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


def parsed_rule_to_rule_dict(parsed: pp.ParseResults):
    name = parsed.name
    forward_rate = parsed.forward_rate
    reverse_rate = parsed.get('reverse_rate', None)

    if not isinstance(name, str):
        raise TypeError("Name must be a string.")

    if not isinstance(forward_rate, str):
        raise TypeError("forward_rate must be a string.")

    if reverse_rate is not None and not isinstance(reverse_rate, str):
        raise TypeError("reverse_rate must be a string or none.")

    left_side = parsed.left_side_reactants
    right_side = parsed.right_side_reactants

    left_reactants = parsed_reactants_to_list(left_side)
    right_reactants = parsed_reactants_to_list(right_side)

    result: ReactionRuleDict = {
        "name": name,
        "reactants": left_reactants,
        "products": right_reactants,
        "forward_rate": forward_rate,
        "reverse_rate": reverse_rate
    }

    return result


def parse_reaction_rule(reaction_str: str):
    """Parses a biochemical reaction and verifies the correctness of complexes and bonds."""

    reagent_pattern = pp.Group(
        COMPLEX_PARSER | pp.Group(MOLECULE_PARSER))
    reagents_pattern = pp.delimitedList(reagent_pattern, '+')

    left_side_reactants = reagents_pattern('left_side_reactants')
    right_side_reactants = reagents_pattern('right_side_reactants')

    arrow_unidirectional = pp.Suppress('->')
    arrow_bidirectional = pp.Suppress('<->')

    expression_parser = pp.Combine(EXPRESSION_PARSER)

    forward_rate_pattern = expression_parser('forward_rate')
    reverse_rate_pattern = expression_parser('reverse_rate')

    name_parser = pp.Optional((NAME_EXPRESSION('name') +
                   pp.Suppress(':')).leaveWhitespace())

    unidirectional_pattern = name_parser + \
        left_side_reactants + \
        arrow_unidirectional + \
        right_side_reactants + \
        forward_rate_pattern

    bidirectional_pattern = name_parser + \
        left_side_reactants + \
        arrow_bidirectional + \
        right_side_reactants + \
        forward_rate_pattern + \
        pp.Literal(',') + \
        reverse_rate_pattern

    rule_parser = unidirectional_pattern | bidirectional_pattern

    try:
        parsed = rule_parser.parseString(reaction_str)
        reaction = parsed_rule_to_rule_dict(parsed)
        return reaction
    except (pp.ParseException, pp.ParseBaseException) as ee:
        raise ParsingError(f"Reaction {reaction_str} not declared correctly.") from ee


def parse_molecule_type(declaration: str):
    """A valid declaration has the format: name() or name(parameters).
    Parameters are comma separated, without spaces after the comma, and can be a word 
    specifying a binding site, or a word followed by one or more tilde and a string 
    specifying a state and possible state values. For example:
        - T(l,Phos~U~P) means T has a binding site with L, and the can be either 
        phosphorylated or unphosphorylated.
        - or T(l,r,Phos~U~P,Meth~A~B~C) same has above, plus a binding side 
        with R and a methilation state with can be either A, B or C
    Returns a dictionary with name, binding_sites, and a dictionary of 
    state dict[str, set[str]], if the name is valid    

    Args:
        declaration (str): _description_

    Returns:
        _type_: _description_
    """

    # Define the format for components enclosed in parentheses
    molecule_type_state = pp.Suppress('~') + STATE
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
        raise ParsingError(
            f"Molecule {declaration} not declared correctly.") from e

    name = parsed.molecule_name
    if not isinstance(name, str):
        raise TypeError("Molecule name must be a string.")

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
    name_pattern = NAME_EXPRESSION('molecule_name')

    pattern_parser = pp.Group(
        COMPLEX_PARSER
        | pp.Group(MOLECULE_PARSER)
        | pp.Group(name_pattern))

    sign_parser = (
        pp.Literal('==')
        | pp.Literal('<=') | pp.Literal('<')
        | pp.Literal('>=') | pp.Literal('>')
    )

    obs_element_parser = pp.Group(
        pattern_parser('pattern')
        + pp.Optional(
            sign_parser('sign') + NUMS('value')
        )
    )

    elements_list_parser = pp.delimitedList(obs_element_parser)

    observable_parser = type_pattern(
        'type') + label_pattern('label') + elements_list_parser('expressions')

    try:
        parsed = observable_parser.parseString(declaration)
        result = parsed_observable_to_dict(parsed)
        return result
    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"Observable '{declaration}' not declared correctly."
        raise ParsingError(msg) from ee


def parse_parameter(declaration: str) -> ParameterDict:
    expression_parser = pp.Combine(EXPRESSION_PARSER)
    parameter_parser = NAME_EXPRESSION('name') + expression_parser('expression')

    try:
        parsed = parameter_parser.parseString(declaration)
        result = parsed_parameter_to_parameter_dict(parsed)
        return result
    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"Parameter {declaration} not declared correctly."
        raise ParsingError(msg) from ee


def parse_seed_species(declaration: str) -> SeedSpeciesDict:
    expression_parser = pp.Combine(EXPRESSION_PARSER)
    reagent_parser = pp.Group(
        COMPLEX_PARSER | pp.Group(MOLECULE_PARSER))

    seed_species_parser = reagent_parser('pattern') + expression_parser('expression')

    try:
        parsed = seed_species_parser.parseString(declaration)
        result = parsed_seed_species_to_seed_species_dict(parsed)
        return result
    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"Seed species {declaration} not declared correctly."
        raise ParsingError(msg) from ee


def parse_compartment(declaration: str) -> CompartmentDict:

    name_parser = NAME_EXPRESSION('name')
    enclosing_compartment_parser = pp.Optional(NAME_EXPRESSION('enclosing_compartment'))
    comment_parser = pp.Optional(COMMENT_PARSER)

    dimensions_parser = pp.Word('23', max=1)('dimensions')
    number_parser = pp.Combine(UNSIGNED_NUMBER_PARSER)

    volume_parser = (number_parser | VARIABLE_PARSER)('volume')

    compartment_parser = (
        name_parser +
        dimensions_parser +
        volume_parser +
        enclosing_compartment_parser +
        comment_parser
    )

    try:
        parsed = compartment_parser.parseString(declaration)
        result = parsed_compartment_to_compartment_dict(parsed)
        return result
    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"Compartment '{declaration}' not declared correctly."
        raise ParsingError(msg) from ee
