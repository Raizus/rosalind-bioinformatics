

from collections import Counter, defaultdict
from itertools import product
from typing import OrderedDict

from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule
from BioInfoToolkit.RuleBasedModel.model.ReactionTransformations import apply_transforms
from BioInfoToolkit.RuleBasedModel.model.Species import Species, species_match_gen
from BioInfoToolkit.RuleBasedModel.network.species_block import SpeciesBlock
from BioInfoToolkit.RuleBasedModel.utils.network_parsers import parse_reaction


class Reaction:
    """Reaction class, to represent a reaction in the network.
    Each element in the list of reactants and products refers to a species id.

    Returns:
        Reaction: _description_
    """
    reactants: list[int]
    products: list[int]
    rate_expression: str
    rule_id: int
    comment: str

    def __init__(self,
                 reactants: list[int],
                 products: list[int],
                 rate_expression: str,
                 rule_id: int,
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

        reaction = Reaction(reactants, products, rate, -1, comment)
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

        match = (set(r_ids) == set(r_ids_other)
                 and set(p_ids) == set(p_ids_other) and
                 self.rate_expression == other.rate_expression)

        return match


def build_rules_dict(rules: OrderedDict[int, ReactionRule]) -> OrderedDict[int, ReactionRule]:
    new_rules: OrderedDict[int, ReactionRule] = OrderedDict()
    count = 1
    for r_id, rule in rules.items():
        if not rule.name:
            rule.name = f"R{r_id}"

        new_rule = rule.get_forward()
        new_rules[count] = new_rule
        count += 1
        if rule.is_bidirectional():
            new_rule = rule.get_reverse()
            new_rules[count] = new_rule
            count += 1
    return new_rules


def breaks_stoich(pattern: Pattern, max_stoich: dict[str, int]):
    mol_counts = pattern.molecule_counts()
    for mol_name, max_count in max_stoich.items():
        count = mol_counts.get(mol_name, 0)
        if count > max_count:
            return True
    return False


def breaks_stoich_multiple(species: list[Pattern], max_stoich: dict[str, int]) -> bool:
    if len(max_stoich) == 0:
        return False

    for specie in species:
        if breaks_stoich(specie, max_stoich):
            return True

    return False


def count_generated_rules(
    react_species_ids: list[int],
    rule_id: int,
    rule: ReactionRule,
    species_block: SpeciesBlock,
    max_stoich: dict[str, int]
) -> defaultdict[tuple[int, tuple[int, ...], tuple[int, ...]], int]:
    """_summary_

    Args:
        react_species_ids (list[int]): _description_
        rule_id (int): _description_
        rule (ReactionRule): _description_
        species_block (SpeciesBlock): _description_
        max_stoich (dict[str, int]): _description_

    Returns:
        defaultdict[tuple[int, tuple[int, ...], tuple[int, ...]], int]: 
        dictionary mapping reactions (tuple of (rule_id, reactants, products)) 
        to counts
    """
    species_dict = species_block.items
    react_species = [
        species_dict[sp_id].pattern for sp_id in react_species_ids]

    reaction_counter: defaultdict[
        tuple[int, tuple[int, ...], tuple[int, ...]],
        int
    ] = defaultdict(int)

    products_gen = apply_transforms(react_species, rule.transformations, rule.modifiers)
    for prod_sp_patts in products_gen:
        # check max stoichiometry
        if breaks_stoich_multiple(prod_sp_patts, max_stoich):
            continue

        # if generated products are new, then add them to species and
        # add reaction to the reaction block
        prod_sp_ids: list[int] = []
        for prod_sp in prod_sp_patts:
            # try to add specie (may not be a new specie)
            specie = Species(prod_sp, "0")
            sp_id = species_block.add_specie(specie)
            prod_sp_ids.append(sp_id)

        reaction_key = (
            rule_id,
            tuple(sorted(react_species_ids)),
            tuple(sorted(prod_sp_ids))
        )
        reaction_counter[reaction_key]+=1

    return reaction_counter


class ReactionGenerator:
    rules: OrderedDict[int, ReactionRule]
    # stores with combinations of rule_ids and sorted reactants have already been computed
    apply_rule_cache: set[tuple[int, tuple[int,...]]]
    # this maps str(pattern): sp_id
    species_match_cache: dict[tuple[str,int], int]

    def __init__(self, rules: OrderedDict[int, ReactionRule]) -> None:
        self.rules = rules
        self.apply_rule_cache = set()
        self.species_match_cache = {}

    def generate(self, species_block: SpeciesBlock, max_stoich: dict[str, int]):
        species_dict = species_block.items

        for rule_id, rule in self.rules.items():
            reactants = rule.reactants

            reactants_gens = [species_match_gen(patt, species_dict,
                                                self.species_match_cache)
                              for patt in reactants]

            reaction_counter: Counter[tuple[int,
                                            tuple[int, ...], tuple[int, ...]]] = Counter()

            for react_sp_ids in product(*reactants_gens):
                rule_sp_key = (rule_id, tuple(sorted(react_sp_ids)))

                # check if we already applied this rule to this combination of reactants
                if rule_sp_key in self.apply_rule_cache:
                    continue

                # apply rule to reactants
                # (id, react_sps, prod_sps)
                counter2 = count_generated_rules(
                    list(react_sp_ids), rule_id, rule,
                    species_block, max_stoich
                )

                reaction_counter.update(counter2)

            for react_key, count in reaction_counter.items():
                _, react_sp_ids, prod_sp_ids = react_key

                rule_sp_key = (rule_id, tuple(sorted(react_sp_ids)))

                # update cache
                self.apply_rule_cache.add(rule_sp_key)

                # create new reaction
                rate_expr = rule.forward_rate
                multiplier = count if not rule.symmetric else 0.5*count
                if multiplier != 1:
                    rate_expr = f"{multiplier}*{rate_expr}"

                comment = f"{rule.name}"
                rxn = Reaction(list(react_sp_ids), list(prod_sp_ids),
                            rate_expr, rule_id, comment)
                yield rxn
