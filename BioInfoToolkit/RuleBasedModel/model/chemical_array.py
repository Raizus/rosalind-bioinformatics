
from functools import partial
from typing import Any
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern


class ChemicalArrayGraph:
    chem_array: list[Pattern]
    graph: nx.Graph

    def __init__(self, chem_array: list[Pattern]) -> None:
        self.chem_array = chem_array


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
        node_map: dict[tuple[int, int], tuple[int, int, int]] = \
            {(j, k): (i, j, k) for j, k in pattern_g.nodes}
        new_pattern_graph: nx.Graph = nx.relabel_nodes(pattern_g, node_map)
        graph: nx.Graph = nx.compose(graph, new_pattern_graph)
    return graph


def node_matching(n1: Any, n2: Any,
                  match_state: bool = False,
                  match_bond: bool = True) -> bool:
    mol1_name = n1.get('molecule_name')
    comp1_name = n1.get('comp_name', None)
    comp1_state = n1.get('state', None)
    comp1_bond = n1.get('bond', None)

    mol2_name = n2.get('molecule_name')
    comp2_name = n2.get('comp_name', None)
    comp2_state = n2.get('state', None)
    comp2_bond = n2.get('bond', None)

    component_match = mol1_name == mol2_name and comp1_name == comp2_name

    # state does not have to match, but must be either undefined in both or defined in both
    state_match = comp1_state == comp2_state
    if not match_state:
        state_match = ((comp1_state == comp2_state)
                       or (comp1_state is not None
                           and comp2_state is not None
                           and len(comp2_state) > 0 and len(comp1_state) > 0))

    # when finding mappings, bond wildcards must match
    is_wildcard = comp1_bond in ('?', '+') or comp2_bond in ('?', '+')
    is_bonded1 = bool(comp1_bond)
    is_bonded2 = bool(comp2_bond)
    bond_match = True
    if is_wildcard:
        bond_match = comp1_bond == comp2_bond
    elif match_bond:
        bond_match = is_bonded1 == is_bonded2

    full_match = component_match and state_match and bond_match
    return full_match


def compare_chemical_array_graphs(graph1: nx.Graph, graph2: nx.Graph) -> bool:
    match_func = partial(node_matching, match_state=True, match_bond=True)
    is_iso: bool = nx.isomorphism.is_isomorphic(graph1, graph2, match_func)
    return is_iso


def compare_chemical_arrays(array1: list[Pattern], array2: list[Pattern]) -> bool:
    graph1 = build_chemical_array_graph(array1)
    graph2 = build_chemical_array_graph(array2)

    return compare_chemical_array_graphs(graph1, graph2)