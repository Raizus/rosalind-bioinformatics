

from itertools import product
from typing import OrderedDict
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule
from BioInfoToolkit.RuleBasedModel.model.Species import Species, species_match_gen
from BioInfoToolkit.RuleBasedModel.network.parsers import parse_reaction


class Reaction:
    """Reaction class, to represent a reaction in the network.
    Each element in the list of reactants and products refers to a species id.

    Returns:
        Reaction: _description_
    """
    reactants: list[int]
    products: list[int]
    rule_id: int
    rate_expression: str
    comment: str

    def __init__(self,
                 reactants: list[int],
                 products: list[int],
                 rule_id: int,
                 rate_expression: str,
                 comment: str = "") -> None:
        """Reaction class for the reaction network

        Args:
            reactants (list[int]): list of reactant species ids
            products (list[int]): list of products species ids
            rule_id (int): rule id that produced this reaction
            rate_expression (str): reaction rate expression
            comment (str, optional): Additional info. Defaults to "".
        """
        self.reactants = reactants
        self.products = products
        self.rule_id = rule_id
        self.rate_expression = rate_expression
        self.comment = comment

    @classmethod
    def from_declaration(cls, declaration: str) -> "Reaction":
        parsed_reaction = parse_reaction(declaration)
        reactants = parsed_reaction['reactants']
        products = parsed_reaction['products']
        rate = parsed_reaction['rate']
        comment = parsed_reaction['comment']

        reaction = Reaction(reactants, products, -1, rate, comment)
        return reaction

    def __repr__(self) -> str:
        reactants_str = ','.join(str(r_id) for r_id in self.reactants)
        products_str = ','.join(str(p_id) for p_id in self.products)

        out = f"{reactants_str} {products_str} {self.rate_expression}"
        if self.comment:
            out += f" # {self.comment}"
        return out

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Reaction):
            return False

        # test equality only by matching species id
        r_ids = self.reactants
        p_ids = self.products

        r_ids_other = other.reactants
        p_ids_other = other.products

        match = (self.rule_id == other.rule_id and
                 set(r_ids) == set(r_ids_other)
                 and set(p_ids) == set(p_ids_other) and
                 self.rate_expression == other.rate_expression)

        return match


def build_rules_dict(rules: OrderedDict[int, ReactionRule]) -> OrderedDict[int, ReactionRule]:
    new_rules: OrderedDict[int, ReactionRule] = OrderedDict()
    count = 1
    for _, rule in rules.items():
        new_rule = rule.get_forward()
        new_rules[count] = new_rule
        count += 1
        if rule.is_bidirectional():
            new_rule = rule.get_reverse()
            new_rules[count] = new_rule
            count += 1
    return new_rules


class ReactionGenerator:
    rules: OrderedDict[int, ReactionRule]
    apply_rule_cache: dict[str, list[int]]
    species_match_cache: dict[tuple[str,int], int]

    def __init__(self, rules: OrderedDict[int, ReactionRule]) -> None:
        self.rules = rules
        self.apply_rule_cache = {}

    def generate(self, species_dict: OrderedDict[int, Species]):
        for rule_id, rule in self.rules.items():
            reactants = rule.reactants
            products = rule.products

            reactants_gens = [species_match_gen(patt, species_dict,
                                                self.species_match_cache)
                              for patt in reactants]

            for react_sp_ids in product(*reactants_gens):
                react_sp_patts = [species_dict[sp_id].pattern for sp_id in react_sp_ids]
                rule_sp_key = (f"({rule_id},"
                               + ','.join(str(sp_id) for sp_id in react_sp_ids)
                               + ')')

                # check if we already applied this rule to this combination of reactants
                if rule_sp_key in self.apply_rule_cache:
                    continue

                # apply rule to reactants
                prod_sp_patts: list[Pattern] = []

                # if generated products are new, then add them to species and 
                # add reaction to the reaction block
                prod_sp_ids: list[int] = []
                for prod_sp in prod_sp_patts:
                    specie = Species(prod_sp, "0")
                    sp_id = species_block.add_species(specie)
                    prod_sp_ids.append(sp_id)

                    # if none of the species is new, we can ignore the reaction
                    comment = f"{rule.name}"
                    rxn = Reaction(list(react_sp_ids), prod_sp_ids, rule_id, comment)
                    # create reaction