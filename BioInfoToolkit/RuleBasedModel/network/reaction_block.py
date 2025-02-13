from BioInfoToolkit.RuleBasedModel.network.blocks import NetworkBlock
from BioInfoToolkit.RuleBasedModel.network.reaction_generation import Reaction
from BioInfoToolkit.RuleBasedModel.utils.utls import format_data_into_lines


from typing import OrderedDict


class ReactionsBlock(NetworkBlock):
    name = "reactions"
    items: OrderedDict[int, Reaction]
    count_id: int

    def __init__(self) -> None:
        self.count_id = 0
        self.items = OrderedDict()

    def add_reaction(self, reaction: Reaction):
        r_id = self.count_id
        self.items[r_id] = reaction
        self.count_id += 1

    def gen_string(self):
        if len(self.items) == 0:
            return ''

        lines: list[str] = []
        lines.append(f"\nbegin {self.name}")

        data: list[tuple[str, str, str, str, str]] = []
        for r_id, reaction in self.items.items():
            comment = f"# {reaction.comment}" if reaction.comment else ''
            reactants_str = (','.join(str(r) for r in reaction.reactants)
                              if len(reaction.reactants) else '0')
            products_str = (','.join(str(r) for r in reaction.products)
                            if len(reaction.products) else '0')
            data_line = (str(r_id),
                         reactants_str,
                         products_str,
                         str(reaction.rate_expression),
                         comment)
            data.append(data_line)

        lines.extend(format_data_into_lines(data))
        lines.append(f"end {self.name}\n")

        return '\n'.join(lines)