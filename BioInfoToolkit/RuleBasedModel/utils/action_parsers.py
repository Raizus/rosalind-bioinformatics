import pyparsing as pp
from typing import Any

from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import NAME_EXPRESSION, UNSIGNED_NUMBER_PARSER, GenerateNetworkDict, ParsingError, SimulateDict




def generate_network_parsed_dict(parsed_dict: dict[Any, Any]) -> GenerateNetworkDict:
    max_iter: int | None = parsed_dict.get('max_iter', None)
    overwrite: bool | None = parsed_dict.get('overwrite', None)
    text_reaction: bool | None = parsed_dict.get('text_reaction', None)
    max_stoich: dict[str, int] | None = parsed_dict.get('max_stoich', None)

    result: GenerateNetworkDict = {
        'max_iter': max_iter,
        'overwrite': overwrite,
        'max_stoich': max_stoich,
        'text_reaction': text_reaction
    }

    return result


def parse_generate_network(declaration: str):
    # Basic patterns
    integer = pp.Word(pp.nums)
    zero_or_one = pp.Word('01', exact=1)
    name = NAME_EXPRESSION

    # parse_actions
    def parse_01_action(token: pp.ParseResults):
        return bool(int(token.asList()[0]))

    def parse_int_action(token: pp.ParseResults):
        return int(token.asList()[0])

    def parse_stoich_action(token: pp.ParseResults):
        stoich_dict: dict[str, int] = {}
        dict_aux = token.asDict()['max_stoich']
        for key, val in dict_aux.items():
            stoich_dict[key] = int(val)
        return stoich_dict

    # Define key=>value pairs for each parameter
    overwrite_expr = (
        pp.Literal("overwrite")
        + pp.Suppress("=>")
        + zero_or_one("overwrite").set_parse_action(parse_01_action))

    text_reaction_expr = (
        pp.Literal("TextReaction")
        + pp.Suppress("=>")
        + zero_or_one("text_reaction").set_parse_action(parse_01_action))

    max_iter_expr = (
        pp.Literal("max_iter")
        + pp.Suppress("=>")
        + integer("max_iter").set_parse_action(parse_int_action))

    # max_stoich list of {name=>value}
    stoich_pair = pp.Group(name("name") + pp.Suppress("=>") +
                        integer("value"))
    max_stoich_expr = (
        pp.Literal("max_stoich") + pp.Suppress("=>")
        + pp.Suppress("{")
        + pp.Dict(pp.delimitedList(stoich_pair)
                  )("max_stoich").set_parse_action(parse_stoich_action)
        + pp.Suppress("}")
    )

    # Combining all expressions
    parameters_expr = pp.Optional(
        pp.Suppress("{") +
        pp.delimitedList(overwrite_expr | text_reaction_expr
                         | max_iter_expr | max_stoich_expr)
        + pp.Suppress("}"))

    expr = (
        pp.Suppress("generate_network(")
        + parameters_expr
        + pp.Suppress(")")
    )

    try:
        parsed = expr.parseString(declaration)
        parsed_dict = parsed.asDict()
        result = generate_network_parsed_dict(parsed_dict)
        return result
    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"generate_network '{declaration}' not declared correctly."
        raise ParsingError(msg) from ee


def parse_simulate(declaration: str) -> SimulateDict:
    # Basic patterns
    integer = pp.Word(pp.nums)
    float_number = pp.Combine(UNSIGNED_NUMBER_PARSER)
    method_string = pp.Suppress('"') + pp.Word(pp.alphas) + pp.Suppress('"')

    # parse_actions
    def parse_method_action(token: pp.ParseResults):
        method = token[0]
        if method != "ssa":
            raise ValueError(
                f"Invalid method '{method}'. Only 'ssa' is supported.")
        return method

    def parse_float_action(token: pp.ParseResults):
        return float(token.asList()[0])

    def parse_int_action(token: pp.ParseResults):
        return int(token.asList()[0])

    # Define key=>value pairs for each parameter
    method_expr = (
        pp.Literal("method")
        + pp.Suppress("=>")
        + method_string("method").set_parse_action(parse_method_action))

    t_end_expr = (
        pp.Literal("t_end")
        + pp.Suppress("=>")
        + float_number("t_end").set_parse_action(parse_float_action))

    t_start_expr = (
        pp.Literal("t_start")
        + pp.Suppress("=>")
        + float_number("t_start").set_parse_action(parse_float_action))

    n_steps_expr = (
        pp.Literal("n_steps")
        + pp.Suppress("=>")
        + integer("n_steps").set_parse_action(parse_int_action))

    # Combining all expressions
    parameters_expr = pp.Optional(
        pp.delimitedList(t_end_expr | t_start_expr | n_steps_expr)
    )

    expr = (
        pp.Suppress("simulate(")
        + pp.Suppress("{")
        + method_expr
        + pp.Suppress(",")
        + parameters_expr
        + pp.Suppress("}")
        + pp.Suppress(")")
    )

    try:
        # Parse the input string
        parsed = expr.parseString(declaration)

        # Convert parsed result to dictionary
        parsed_dict = parsed.asDict()

        # Set defaults for optional parameters
        t_start = parsed_dict.get("t_start", 0.0)  # default is 0
        n_steps = parsed_dict.get("n_steps", None)  # default is None

        # Ensure required parameters are present
        t_end = parsed_dict["t_end"]  # t_end is mandatory

        # Build the result dictionary
        result: SimulateDict = {
            "method": parsed_dict["method"],
            "t_start": t_start,
            "t_end": t_end,
            "n_steps": n_steps
        }

        return result

    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"simulate '{declaration}' not declared correctly."
        raise ParsingError(msg) from ee