from typing import Any
import pyparsing as pp

from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import COMMENT_PARSER, \
    COMPLEX_PARSER, EXPRESSION_PARSER, LABEL_PARSER, MOLECULE_PARSER, \
    NUMS, VARIABLE_PARSER, NetworkParameterDict, NetworkSeedSpeciesDict, ObservablesGroupDict, ParsingError, \
    ReactionDict, parsed_parameter_to_parameter_dict, parsed_seed_species_to_seed_species_dict


def parsed_obs_group_to_dict(parsed: pp.ParseResults | Any) -> ObservablesGroupDict:
    if not isinstance(parsed, pp.ParseResults):
        raise TypeError(f"{parsed} must be of type ParseResults.")

    g_id = parsed.g_id
    if not isinstance(g_id, str):
        raise TypeError(f"type {g_id} must be a string")

    label = parsed.label
    if not isinstance(label, str):
        raise TypeError(f"type {label} must be a string")

    parsed_species = parsed.species
    species: list[tuple[int, int]] = []
    for parsed_specie in parsed_species:
        if not isinstance(parsed_specie, pp.ParseResults):
            raise TypeError(f"{parsed_specie} must be of type ParseResults.")

        sp_id = parsed_specie.sp_id
        if not isinstance(sp_id, str):
            raise TypeError(f"type {sp_id} must be a string")
        
        weight = parsed_specie.get('weight', "1")
        if not isinstance(weight, str):
            raise TypeError(f"type {weight} must be a string")

        species.append((int(sp_id), int(weight)))

    result: ObservablesGroupDict = {
        "id": int(g_id),
        "label": label,
        "species": species
    }
    return result


def parsed_reaction_to_dict(parsed: pp.ParseResults | Any) -> ReactionDict:
    if not isinstance(parsed, pp.ParseResults):
        raise TypeError(f"{parsed} must be of type ParseResults.")

    r_id = parsed.r_id
    if not isinstance(r_id, str):
        raise TypeError(f"type {r_id} must be a string")

    parsed_reactants = parsed.reactants
    if not isinstance(parsed_reactants, pp.ParseResults):
        raise TypeError(
            f"type {parsed_reactants} must be of type ParseResults")

    reactants: list[int] = []
    for val in parsed_reactants.as_list():
        if int(val) == 0:
            continue
        reactants.append(int(val))

    parsed_products = parsed.products
    if not isinstance(parsed_products, pp.ParseResults):
        raise TypeError(f"type {parsed_products} must be of type ParseResults")

    products: list[int] = []
    for val in parsed_products.as_list():
        if int(val) == 0:
            continue
        products.append(int(val))

    rate = parsed.rate
    if not isinstance(rate, str):
        raise TypeError(f"type {rate} must be a string")

    comment = parsed.get('comment', '')
    if not isinstance(comment, str):
        raise TypeError(f"type {comment} must be a string")

    result: ReactionDict = {
        "id": int(r_id),
        "reactants": reactants,
        "products": products,
        "rate": rate,
        "comment": comment
    }
    return result


def parse_reaction(declaration: str) -> ReactionDict:
    """
    "{g_id} {name}"

    Args:
        declaration (str): _description_
    """

    id_parser = NUMS('r_id')
    species_list_parser = pp.delimited_list(NUMS)
    reactants_parser = species_list_parser('reactants')
    products_parser = species_list_parser('products')
    expression_parser = pp.Combine(EXPRESSION_PARSER)
    comment_parser = pp.Optional(COMMENT_PARSER)

    reaction_parser = (id_parser + reactants_parser + products_parser +
                       expression_parser('rate') + comment_parser)

    try:
        parsed = reaction_parser.parseString(declaration)
    except (pp.ParseException, pp.ParseBaseException) as exc:
        msg = f"Reaction declaration '{declaration}' is not in the correct format."
        raise ParsingError(msg) from exc

    result = parsed_reaction_to_dict(parsed)
    return result


def parse_observables_group(declaration: str) -> ObservablesGroupDict:
    """
    "{g_id} {name}"

    Args:
        declaration (str): _description_
    """

    id_parser = NUMS('g_id')
    name_parser = LABEL_PARSER('label')
    element_parser = pp.Optional(
        NUMS('weight') + pp.Suppress('*')) + NUMS('sp_id')
    elements_parser = pp.delimited_list(pp.Group(element_parser))

    obs_group_parser = id_parser + name_parser + pp.Optional(elements_parser('species'))

    try:
        parsed = obs_group_parser.parseString(declaration)
    except (pp.ParseException, pp.ParseBaseException) as exc:
        msg = f"Observables group declaration '{declaration}' is not in the correct format."
        raise ParsingError(msg) from exc

    result = parsed_obs_group_to_dict(parsed)
    return result


def parse_parameters(declaration: str) -> NetworkParameterDict:
    id_parser = NUMS('id')
    expression_parser = pp.Combine(EXPRESSION_PARSER)
    comment_parser = pp.Optional(COMMENT_PARSER)
    parameter_parser = (id_parser +
                        VARIABLE_PARSER('name') +
                        expression_parser('expression') +
                        comment_parser)

    try:
        parsed = parameter_parser.parseString(declaration)

    except (pp.ParseException, pp.ParseBaseException) as ee:
        raise ParsingError(
            f"Parameter {declaration} not declared correctly.") from ee

    p_id = parsed.id
    if not isinstance(p_id, str):
        raise ValueError("Parameter p_id must be a string.")

    parsed_param = parsed_parameter_to_parameter_dict(parsed)

    result: NetworkParameterDict = {
        'id': int(p_id),
        **parsed_param
    }
    return result


def parse_seed_species(declaration: str) -> NetworkSeedSpeciesDict:
    expression_parser = pp.Combine(EXPRESSION_PARSER)
    reagent_parser = pp.Group(
        COMPLEX_PARSER | pp.Group(MOLECULE_PARSER))
    id_parser = NUMS('sp_id')

    seed_species_parser = (id_parser +
                           reagent_parser('pattern') +
                           expression_parser('expression'))

    try:
        parsed = seed_species_parser.parseString(declaration)
    except (pp.ParseException, pp.ParseBaseException) as ee:
        raise ParsingError(
            f"Seed species {declaration} not declared correctly.") from ee

    sp_id = parsed.sp_id
    if not isinstance(sp_id, str):
        raise ValueError("Parameter sp_id must be a string.")

    parsed_species = parsed_seed_species_to_seed_species_dict(parsed)

    result: NetworkSeedSpeciesDict = {
        'id': int(sp_id),
        **parsed_species
    }
    return result
