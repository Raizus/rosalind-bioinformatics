from collections import defaultdict
from functools import partial
from typing import Any
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Parsers import parse_reaction
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern
from BioInfoToolkit.RuleBasedModel.model.ReactionTransformations import ChangeStateAction, ReactionTransformation, apply_transforms


def node_matching(n1: Any, n2: Any, 
                  match_state: bool = False,
                  match_bond: bool = True) -> bool:
    label1: tuple[str, str | None, str |
                  None, str | None] = n1['label']
    label2: tuple[str, str | None, str |
                  None, str | None] = n2['label']
    mol1_name, comp1_name, comp1_state, comp1_bond = label1
    mol2_name, comp2_name, comp2_state, comp2_bond = label2

    component_match = mol1_name == mol2_name and comp1_name == comp2_name
    # state does not have to match, but must be either undefined in both or defined in both
    state_match = comp1_state == comp2_state
    if not match_state:
        state_match = ((comp1_state == comp2_state)
                       or (comp1_state is not None
                           and comp2_state is not None
                           and len(comp2_state) > 0 and len(comp1_state) > 0))

    # check bonds
    comp1_is_bonded = (comp1_bond is not None
                       and len(comp1_bond) > 0
                       and comp1_bond != '?')
    comp2_is_bonded = (comp2_bond is not None
                       and len(comp2_bond) > 0
                       and comp2_bond != '?')

    bond_match = True
    if match_bond:
        bond_match = comp1_is_bonded == comp2_is_bonded
        if comp1_bond == '?' or comp2_bond == '?':
            bond_match = True
        elif comp1_bond == '+':
            bond_match = comp2_bond is not None and len(comp2_bond) > 0

    full_match = component_match and state_match and bond_match
    return full_match


def build_chemical_array_graph(chemical_array: list[Pattern]) -> nx.Graph:
    """If we build the reactants graph (which is equivalent to combining the graphs of each pattern,
    with nodes relabeled so patterns don't combine), 
    then we can analize the graphs of reactants and products

    Args:
        reactants (list[Pattern]): _description_
    """
    graph = nx.Graph()
    # create graph
    for i, pattern in enumerate(chemical_array):
        pattern_g = pattern.graph
        node_map: dict[tuple[int,int], tuple[int,int,int]] = \
            {(j,k): (i,j,k) for j,k in pattern_g.nodes}
        new_pattern_graph: nx.Graph = nx.relabel_nodes(pattern_g, node_map)
        graph: nx.Graph = nx.compose(graph, new_pattern_graph)
    return graph


def compare_chemical_array_graphs(graph1: nx.Graph, graph2: nx.Graph) -> bool:
    match_func = partial(node_matching, match_state=True)
    is_iso: bool = nx.isomorphism.is_isomorphic(graph1, graph2, match_func)
    return is_iso


class ReactionRule:
    """Allowed transformations:
        - Forming a bond:  A(b) + B(a) -> A(b!0).B(a!0)
        - Breaking a bond:  A(b!0).B(a!0) -> A(b) + B(a)
        - Changing of component state:  X(y\~0) -> X(y\~p)
        - Creating a molecule:  A(b) -> A(b) + C(d)
        - Destroying a molecule:  A(b) + B(a) -> A(b)

    A single reaction may involve any number of transformations:
        - A(b) + B(a!0).C(d!0) -> A(b!0).B(a!0) + C(d)  Break bond and form new one
        - A(b) + B(a!0).C(d!0) -> A(b!0).B(a!0)         Break bond, form new bond, destroy C molecule

    A Reaction rule is a composition of simple transformations, that transforms the reactants 
    patterns into the product patterns (and vice-versa if the reaction is bidirectional)

    If we create a graph of the reactants and products, then the reaction rule 
    can be seen as a composition of transformations that transforms the reactants graph 
    into the products graph.
        - Forming a bond -> add an edge between two molecule component nodes, 
        and change the bond label for the two involved components
        - Removing a bond -> removing an edge between two molecule component nodes, 
        and change the bond label for the two involved components
        - Change component state -> change the label of the component state
        - Creating a molecule -> Add a molecule node and adjacent component nodes 
        (connect molecule to each component)
        - Destroying a molecule -> Remove a molecule node and adjacent component nodes 
        (connect molecule to each component). Can only destroy a molecule if it's graph 
        (molecule + component nodes) is not connected to any other nodes, 
        meaning if there is a bond, first we need to break the bond, 
        and only after can we destroy a molecule 
    """

    name: str
    reactants: list[Pattern] = []
    products: list[Pattern] = []
    forward_rate: str # forward reaction rate expression
    forward_rate_value: float
    reverse_rate: str|None  # forward reaction rate expression
    reverse_rate_value: float
    reactants_graph: nx.Graph
    products_graph: nx.Graph

    def __init__(self,
                 name: str,
                 reactants: list[Pattern],
                 products: list[Pattern],
                 forward_rate: str,
                 reverse_rate: str|None = None) -> None:
        self.name = name
        self.forward_rate = forward_rate
        self.reactants = reactants
        self.products = products
        self.forward_rate_value = 0.0
        self.reverse_rate = reverse_rate
        self.reverse_rate_value = 0.0

        self.reactants_graph = build_chemical_array_graph(reactants)
        self.products_graph = build_chemical_array_graph(products)

    @classmethod
    def from_declaration(cls, declaration: str, molecules: dict[str, MoleculeType]):
        reaction_dict = parse_reaction(declaration)

        if reaction_dict is None:
            raise ValueError(f"Invalid reaction declaration: {declaration}")

        name = reaction_dict["name"]
        forward_rate = reaction_dict["forward_rate"]
        reverse_rate = reaction_dict["reverse_rate"]
        parsed_reactants = reaction_dict["reactants"]
        parsed_products = reaction_dict["products"]
        reactants = [Pattern.from_dict(parsed_pattern, molecules) 
                     for parsed_pattern in parsed_reactants]
        products = [Pattern.from_dict(parsed_pattern, molecules)
                     for parsed_pattern in parsed_products]
        reaction = ReactionRule(name, reactants, products, forward_rate, reverse_rate)

        return reaction

    def is_bidirectional(self) -> bool:
        return self.reverse_rate is not None

    def __repr__(self) -> str:
        left_side = ' + '.join(str(reagent) for reagent in self.reactants)
        right_side = ' + '.join(str(reagent) for reagent in self.products)
        arrow = '<->' if self.is_bidirectional() else '->'

        rates_str = f"{self.forward_rate}"
        if self.is_bidirectional():
            rates_str += ', ' + f"{self.reverse_rate}"

        out = f"{self.name}: {left_side} {arrow} {right_side} {rates_str}"
        return out


def reactants_tuples_gen(reactants: list[Pattern]):
    for i, pattern in enumerate(reactants):
        for j, molecule in enumerate(pattern.molecules):
            yield ((i,j), (pattern, molecule))


def _find_pattern_to_pattern_mapping(
        reactants: list[Pattern],
        products: list[Pattern]
) -> tuple[dict[int, set[int]], dict[int, set[int]]]:
    """ Finds a map between reagents and products
    A reactant matches with a product if their graphs are isomorphic (ignoring state changes)
    This can be used to detect state changes in components.

    Args:
        reactants (list[Pattern]): _description_
        products (list[Pattern]): _description_

    Returns:
        dict[int, set[int]]: _description_
    """
    # finds a map between reagents and products
    # a reactant matches with a product if their graphs are isomorphic (ignoring state changes)
    # this can be used to detect state changes in components.

    # reactants_g = build_reactants_graph(reactants)
    # products_g = build_reactants_graph(products)

    mapping: dict[int, set[int]] = {}
    for i, reactant in enumerate(reactants):
        mapping[i] = set()
        for j, product in enumerate(products):
            g1 = reactant.graph
            g2 = product.graph
            is_iso: bool = nx.isomorphism.is_isomorphic(g1, g2, node_matching)
            if is_iso:
                mapping[i].add(j)

    reverse_map: dict[int, set[int]] = {}
    for i, js in mapping.items():
        for j in js:
            if j in reverse_map:
                reverse_map[j].add(i)
            else:
                reverse_map[j] = set([i])

    return mapping, reverse_map


def find_molecule_to_molecule_mapping(reactants: list[Pattern],
                                      products: list[Pattern]):

    for n1_id, (patt1, mol1) in reactants_tuples_gen(reactants):
        for n2_id, (patt2, mol2) in reactants_tuples_gen(products):
            # check if molecule graphs are isomorphic
            match_func = partial(node_matching, match_state=False)
            nx.isomorphism.is_isomorphic(g1, g2, )
            pass


def _find_molecule_mapping(
        reactants: list[Pattern],
        products: list[Pattern]
        ) -> dict[tuple[int, int], set[tuple[int, int]]]:
    # finds a map between reagents and products
    # the keys of the dictionary are a tuple (i,j) means the molecule j of reactant i,
    # the same applies for values (molecule j of product i)

    # reactants -> products
    mapping_in_out: dict[tuple[int, int], set[tuple[int, int]]] = {}
    for (n1_id, (_, molecule1)) in reactants_tuples_gen(reactants):
        # find matching molecules in the products
        mapping_in_out[n1_id] = set()
        comp_counts1 = molecule1.components_counts
        for (n2_id, (_, molecule2)) in reactants_tuples_gen(products):
            comp_counts2 = molecule2.components_counts
            if molecule1.name == molecule2.name and comp_counts1 == comp_counts2:
                mapping_in_out[n1_id].add(n2_id)

    return mapping_in_out


def find_molecule_mappings(reactants: list[Pattern], products: list[Pattern]):
    aux_mapping_in_out = _find_molecule_mapping(reactants, products)
    aux_mapping_out_in = _find_molecule_mapping(reactants=products, products=reactants)

    # we want to build a map where each molecule maps to at most 1 other molecule
    mapping_in_out: dict[tuple[int, int], tuple[int, int]] = {}
    mapping_out_in: dict[tuple[int, int], tuple[int, int]] = {}

    # We want each molecule in the reactants to map to at most one other
    # molecule in the products.
    # However, if a molecule is deleted then a molecule in the reactants side
    # does not map to any molecule in the products side.

    # assign trivial mappings first, then the rest

    for n_id_in, y in aux_mapping_in_out.items():
        if len(y) == 1:
            n_id_out = list(y)[0]
            y2 = mapping_out_in.get(n_id_out, {})
            if len(y2) == 1:
                n_id_out2 = list(y2)[0]
                if (n_id_out2 == n_id_in):
                    mapping_in_out[n_id_in] = n_id_out
                    mapping_out_in[n_id_in] = n_id_out

    # update provisory mappings
    for n1 in mapping_in_out:
        aux_mapping_in_out.pop(n1)
    for n1 in mapping_out_in:
        aux_mapping_out_in.pop(n1)

    # now the other cases
    # if there's a value in aux_mapping_in_out with length of 0, then molecules are deleted
    # if there's a value in aux_mapping_out_in with length of 0, then molecules are created

    # we are left with other mappings


def find_bijective_map(mapping: dict[int, set[int]], reverse_mapping: dict[int, set[int]]):
    bijective_map: dict[int, int] = {}

    for i, js in mapping.items():
        if len(js) == 1:
            j = next(iter(js))
            _is = reverse_mapping[j]
            if len(_is) == 1:
                i2 = next(iter(js))
                if i2 == i:
                    bijective_map[i] = j
    return bijective_map


def find_state_changes(reactants: list[Pattern], products: list[Pattern]):
    reactants_map, reverse_map = _find_pattern_to_pattern_mapping(
        reactants, products)
    bijective_map = find_bijective_map(reactants_map, reverse_map)

    considered_reactants: set[int] = set()
    considered_products: set[int] = set()

    transformations: list[ReactionTransformation] = []

    for i, js in reactants_map.items():
        # len(js) == 0 this can mean created / destroyed molecules
        # or bond forming / breaking, ignore for now

        # we have a 1 to 1 reactant to product map (check for state changes)
        if len(js) == 1 and i in bijective_map and bijective_map[i] == next(iter(js)):
            j = next(iter(js))
            reactant = reactants[i]
            product = products[j]

            # either state changes or no changes at all
            # are isomorphic -> no changes to this reactant from input to output
            if reactant == product:
                considered_reactants.add(i)
                considered_products.add(j)
            else:  # 1 or more state changes
                # find component changes
                # check all isomorphism mappings, and compare component by component
                matcher = nx.isomorphism.GraphMatcher(
                    reactant.graph, product.graph, node_matching)
                for node_map in matcher.isomorphisms_iter():
                    for n1, n2 in node_map.items():
                        _, comp1_name, comp1_state, _ = reactant.graph.nodes[n1]['label']
                        _, comp2_name, comp2_state, _ = product.graph.nodes[n2]['label']

                        if comp1_name and comp2_name and comp1_state != comp2_state:
                            center_in: list[tuple[int, int, int]] = [(i, *n1)]
                            center_out: list[tuple[int, int, int]] = [(j, *n2)]
                            action = ChangeStateAction(
                                reactants, products, center_in, center_out)
                            transformations.append(action)
                    break

    return transformations, considered_reactants, considered_products


def find_molecule_creation(reactants: list[Pattern], products: list[Pattern]):
    pass


def decompose_reaction(reaction: ReactionRule) -> list[ReactionTransformation]:
    reactants = reaction.reactants
    products = reaction.products
    current_reactants = reactants
    all_transformations: list[ReactionTransformation] = []

    while not compare_chemical_array_graphs(
        build_chemical_array_graph(current_reactants),
                                   reaction.products_graph):
        # state change transformations do not change the number or order of reactants
        transformations, considered_reactants, considered_products = \
            find_state_changes(current_reactants, products)
        
        if len(transformations):
            all_transformations.extend(transformations)
            current_reactants = apply_transforms(current_reactants, transformations)
            continue
        
        # apply transformations and check new products, if it's the same we're finished

        # Note that creating and destroying molecules or breaking and forming bonds will change the number of reactants in a list

        # Start with molecule creation and destruction

        # Then bond breaking

        # Then bond forming

    return transformations


def verify_transformations(reactants: list[Pattern],
                           transformations: list[ReactionTransformation],
                           products_g: nx.Graph):
    """Checks if the transformations are applied to reactants, 
    the resulting products graph matches the given products_g

    Args:
        reactants (list[Pattern]): _description_
        transformations (list[ReactionTransformation]): _description_
        products_g (nx.Graph): _description_

    Returns:
        _type_: _description_
    """
    products_res = apply_transforms(reactants, transformations)
    prods2_g = build_chemical_array_graph(products_res)
    g_equal = compare_chemical_array_graphs(products_g, prods2_g)
    return g_equal
