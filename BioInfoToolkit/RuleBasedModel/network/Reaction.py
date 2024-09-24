from BioInfoToolkit.RuleBasedModel.model.Species import Species


class Reaction:
    reactants: list[Species]
    products: list[Species]
    reaction_id: int
    rule_id: int
    rate_expression: str
    comment: str

    def __init__(self, reactants: list[Species],
                 products: list[Species],
                 reaction_id: int,
                 rule_id: int) -> None:
        self.reactants = reactants
        self.products = products
        self.reaction_id = reaction_id
        self.rule_id = rule_id
        self.comment = ""

    def __repr__(self) -> str:
        reactants_str = ','.join(str(reactant.id) for reactant in self.reactants)
        products_str = ','.join(str(product.id)
                                 for product in self.products)

        out = f"{self.reaction_id} {reactants_str} {products_str} {self.rate_expression}"
        if self.comment:
            out += f" {self.comment}"
        return out

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Reaction):
            return False

        # test equality only by matching species id
        r_ids = [reactant.id for reactant in self.reactants]
        p_ids = [product.id for product in self.products]

        r_ids_other = [reactant.id for reactant in other.reactants]
        p_ids_other = [product.id for product in other.products]

        match = (set(r_ids) == set(r_ids_other)
                 and set(p_ids) == set(p_ids_other) and
                 self.rate_expression == other.rate_expression)

        return match