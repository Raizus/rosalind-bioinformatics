import pyparsing as pp
from typing import Any, TypedDict

from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import COMPLEX_PARSER, EXPRESSION_PARSER, MOLECULE_PARSER, NAME_EXPRESSION, UNSIGNED_NUMBER_PARSER, MoleculeDict, ParsingError, parsed_pattern_to_dict_list
from BioInfoToolkit.RuleBasedModel.utils.utls import eval_expr


class GenerateNetworkDict(TypedDict):
    overwrite: bool | None
    text_reaction: bool | None
    max_stoich: dict[str, int] | None
    max_iter: int | None


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


class SimulateDict(TypedDict):
    method: str
    t_start: float
    t_end: float
    n_steps: int | None
    continue_: bool | None
    verbose: bool | None
    atol: float | None
    rtol: float | None
    tau: float | None


def action_float_param_parser(param_name: str, parsed_name: str):
    float_number = pp.Combine(UNSIGNED_NUMBER_PARSER)

    def parse_float_param(token: pp.ParseResults):
        return float(token.asList()[0])

    expr = (
        pp.Literal(param_name)
        + pp.Suppress("=>")
        + float_number(parsed_name).set_parse_action(parse_float_param))

    return expr


def action_01_param_parser(param_name: str, parsed_name: str):
    zero_or_one = pp.Word('01', exact=1)

    def parse_01_param(token: pp.ParseResults):
        return bool(int(token.asList()[0]))

    expr = (
        pp.Literal(param_name)
        + pp.Suppress("=>")
        + zero_or_one(parsed_name).set_parse_action(parse_01_param))

    return expr


def parse_simulate(declaration: str) -> SimulateDict:
    # Basic patterns
    integer = pp.Word(pp.nums)
    exp_part = pp.Optional(pp.CaselessLiteral(
        'e') + pp.Word(pp.nums))
    num_expr = pp.Combine(integer + exp_part)
    method_string = pp.Suppress('"') + pp.Word(pp.alphas, pp.alphas+'-') + pp.Suppress('"')

    # parse_actions
    def parse_method_action(token: pp.ParseResults):
        method = token[0]
        if method not in ("ssa", "ode", "tau-leap"):
            raise ValueError(
                f"Invalid method '{method}'. Only 'ssa', 'ode' and 'tau-leap' is supported.")
        return method

    def parse_int_action(token: pp.ParseResults):
        string: str = token.asList()[0]
        _dict: dict[str, int|float] = {}
        res, _ = eval_expr(string, _dict)
        if isinstance(res, (int, float)):
            return int(res)
        raise ValueError(f"Parameter with value '{string}' must evaluate to integer.")

    # Define key=>value pairs for each parameter
    method_expr = (
        pp.Literal("method")
        + pp.Suppress("=>")
        + method_string("method").set_parse_action(parse_method_action))

    t_end_expr = action_float_param_parser("t_end", "t_end")
    t_start_expr = action_float_param_parser("t_start", "t_start")
    atol_expr = action_float_param_parser("atol", "atol")
    rtol_expr = action_float_param_parser("rtol", "rtol")
    tau_expr = action_float_param_parser("tau", "tau")
    
    continue_expr = action_01_param_parser("continue", "continue_")
    verbose_expr = action_01_param_parser("verbose", "verbose")

    n_steps_expr = (
        pp.Literal("n_steps")
        + pp.Suppress("=>")
        + num_expr("n_steps").set_parse_action(parse_int_action))

    # Combining all expressions
    parameters_expr = pp.Optional(
        pp.delimitedList(t_end_expr | t_start_expr | n_steps_expr
                         | continue_expr | atol_expr | rtol_expr
                         | tau_expr | verbose_expr)
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
        continue_val = parsed_dict.get("continue_", None)  # default is None
        atol = parsed_dict.get("atol", None)
        rtol = parsed_dict.get("rtol", None)
        tau = parsed_dict.get("tau", None)
        verbose = parsed_dict.get("verbose", None)

        # Ensure required parameters are present
        t_end = parsed_dict["t_end"]  # t_end is mandatory

        # Build the result dictionary
        result: SimulateDict = {
            "method": parsed_dict["method"],
            "t_start": t_start,
            "t_end": t_end,
            "n_steps": n_steps,
            "continue_": continue_val,
            "atol": atol,
            "rtol": rtol,
            "tau": tau,
            "verbose": verbose,
        }

        return result

    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"simulate '{declaration}' not declared correctly."
        raise ParsingError(msg) from ee


class SetConcentrationDict(TypedDict):
    pattern: list[MoleculeDict]
    expression: str


def parsed_set_concentration_to_dict(parsed: pp.ParseResults) -> SetConcentrationDict:
    parsed_pattern = parsed.species
    pattern = parsed_pattern_to_dict_list(parsed_pattern)

    expression = parsed.expression
    if not isinstance(expression, str):
        raise ValueError("See species expression must be a string.")

    result: SetConcentrationDict = {
        'pattern': pattern,
        'expression': expression,
    }
    return result


def parse_set_concentration(declaration: str) -> SetConcentrationDict:
    # Basic patterns
    expression_parser = pp.Combine(EXPRESSION_PARSER)
    reagent_parser = pp.Group(
        COMPLEX_PARSER | pp.Group(MOLECULE_PARSER))

    # Define key=>value pairs for each parameter
    species_parser = (
        pp.Suppress('"') + reagent_parser('species') + pp.Suppress('"')
    )

    concentration_parser = (
        pp.Suppress('"') + expression_parser('expression') + pp.Suppress('"')
    )

    # Combining all expressions
    set_conc_parser = (
        pp.Suppress("setConcentration(")
        + species_parser
        + pp.Suppress(',')
        + concentration_parser
        + pp.Suppress(')')
    )

    try:
        parsed = set_conc_parser.parseString(declaration)
        result = parsed_set_concentration_to_dict(parsed)
        return result
    except (pp.ParseException, pp.ParseBaseException) as ee:
        msg = f"setConcentration '{declaration}' not declared correctly."
        raise ParsingError(msg) from ee
