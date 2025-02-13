from collections import defaultdict
from functools import partial
from typing import Any, Generator
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.compartment import Compartments
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.RuleModifiers import RuleModifiers
from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_reaction_rule
from BioInfoToolkit.RuleBasedModel.model.Pattern import Molecule, Pattern
from BioInfoToolkit.RuleBasedModel.model.ReactionTransformations import BreakBondTransform, \
    ChangeStateAction, CreateMoleculeAction, DestroyMoleculeAction, FormBondTransform, \
    ReactionTransformation, apply_transforms
from BioInfoToolkit.RuleBasedModel.model.chemical_array import build_chemical_array_graph, \
    compare_chemical_array_graphs, node_matching
from BioInfoToolkit.RuleBasedModel.utils.utls import eval_expr

ChemArrayNodeId = tuple[int, int, int]

def chem_array_graph_automorphism_gen(graph: nx.Graph):
    match_func = partial(node_matching, match_state=True, match_bond=True)
    matcher = nx.isomorphism.GraphMatcher(graph, graph, match_func)
    yield from matcher.isomorphisms_iter()


class ReactionRule:
    """Allowed transformations:
        - Forming a bond:  A(b) + B(a) -> A(b!0).B(a!0)
        - Breaking a bond:  A(b!0).B(a!0) -> A(b) + B(a)
        - Changing of component state:  X(y\~0) -> X(y\~p)
        - Creating a molecule:  A(b) -> A(b) + C(d)
        - Destroying a molecule:  A(b) + B(a) -> A(b)

    A single reaction may involve any number of transformations:
        - A(b) + B(a!0).C(d!0) -> A(b!0).B(a!0) + C(d)  Break bond and form new one
        - A(b) + B(a!0).C(d!0) -> A(b!0).B(a!0)         Break bond, form new bond, 
                                                        destroy C molecule

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
    transformations: list[ReactionTransformation]
    modifiers: RuleModifiers | None
    symmetric: bool

    def __init__(self,
                 name: str,
                 reactants: list[Pattern],
                 products: list[Pattern],
                 forward_rate: str,
                 reverse_rate: str|None = None,
                 modifiers: RuleModifiers | None = None) -> None:
        self.name = name
        self.forward_rate = forward_rate
        self.reactants = reactants
        self.products = products
        self.forward_rate_value = 0.0
        self.reverse_rate = reverse_rate
        self.reverse_rate_value = 0.0
        self.transformations = []
        self.modifiers = modifiers

        self.reactants_graph = build_chemical_array_graph(reactants)
        self.products_graph = build_chemical_array_graph(products)
        self.symmetric = self.is_symmetric()

    @classmethod
    def from_declaration(cls, declaration: str, molecules: dict[str, MoleculeType]):
        reaction_dict = parse_reaction_rule(declaration)

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
        parsed_mods = reaction_dict['modifiers']
        modifiers = RuleModifiers(parsed_mods)
        reaction = ReactionRule(name, reactants, products, forward_rate, reverse_rate, modifiers)

        return reaction

    def is_symmetric(self) -> bool:
        def automorphism_counter(iso_gen: Generator[Any, Any, None]):
            count = 0
            # we only want to count automorphisms where molecules map to different molecules.
            # Consider a molecule with two components of the same type. This molecule will
            # have two automorphisms but only one should count to the determine
            # if the reaction is symmetric
            mol_maps: list[dict[tuple[int, int], tuple[int, int]]] = []
            for node_map in iso_gen:
                mol_map: dict[tuple[int, int], tuple[int, int]] = \
                    {n1[:2]: n2[:2] for n1, n2 in node_map.items()
                           if n1[2] == -1 and n2[2] == -1}
                if mol_map not in mol_maps:
                    mol_maps.append(mol_map)
                    count += 1
            return count

        gen1 = chem_array_graph_automorphism_gen(self.reactants_graph)
        react_symmetry_c = automorphism_counter(gen1)
        gen2 = chem_array_graph_automorphism_gen(self.products_graph)
        prod_symmetry_c = automorphism_counter(gen2)
        if (react_symmetry_c > 1 and prod_symmetry_c > 1
            and react_symmetry_c == prod_symmetry_c):
            return True
        return False

    def is_bidirectional(self) -> bool:
        return self.reverse_rate is not None

    def get_forward(self) -> "ReactionRule":
        new_rule = ReactionRule(self.name, self.reactants, self.products,
                                self.forward_rate, None, self.modifiers)
        new_rule.transformations = self.transformations
        return new_rule

    def get_reverse(self) -> "ReactionRule":
        if not self.reverse_rate:
            raise TypeError(f"Reaction rule must have a reverse rate ({self.reverse_rate})")

        name = f"_reverse_{self.name}"
        new_rule = ReactionRule(name, self.products, self.reactants,
                                self.reverse_rate, None, self.modifiers)
        new_rule.decompose_reaction()
        return new_rule

    def validate_reactants(
        self,
        molecule_types: dict[str, MoleculeType],
        compartments: Compartments
    ) -> bool:
        for pattern in self.reactants + self.products:
            valid = pattern.validate(molecule_types, compartments)
            if not valid:
                return False

        return True

    def validate_rates(self, variables: dict[str, int|float]) -> bool:
        rate_f = self.forward_rate
        rate_r = self.reverse_rate

        value, _ = eval_expr(rate_f, variables)
        if not isinstance(value, int) and not isinstance(value, float):
            return False

        if rate_r:
            value, _ = eval_expr(rate_r, variables)
            if not isinstance(value, int) and not isinstance(value, float):
                return False

        return True

    def decompose_reaction(self):
        transformations = decompose_reaction(self)
        self.transformations = transformations

    def __repr__(self) -> str:
        left_side = ' + '.join(str(reagent) for reagent in self.reactants) \
            if len(self.reactants) > 0 else '0'
        right_side = ' + '.join(str(reagent) for reagent in self.products) \
            if len(self.products) > 0 else '0'
        arrow = '<->' if self.is_bidirectional() else '->'

        rates_str = f"{self.forward_rate}"
        if self.is_bidirectional():
            rates_str += ', ' + f"{self.reverse_rate}"

        out = f"{left_side} {arrow} {right_side} {rates_str}"
        if self.name:
            out = f"{self.name}: {out}"
        return out


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
    products_res = next(apply_transforms(reactants, transformations))
    prods2_g = build_chemical_array_graph(products_res)
    g_equal = compare_chemical_array_graphs(products_g, prods2_g)
    return g_equal


def reactant_idxs_iter(
    chemical_array: list[Pattern]
) -> Generator[list[ChemArrayNodeId], Any, None]:
    """Iterator that yields a list of node ids for each molecule in the chemical array

    Args:
        chemical_array (list[Pattern]): _description_

    Yields:
        Generator[list[ChemArrayNodeId], Any, None]: list of node ids for each molecule
    """
    for i, reactant in enumerate(chemical_array):
        for j, mol in enumerate(reactant.molecules):
            idxs: list[ChemArrayNodeId] = [(i, j, -1)]
            for k, _ in enumerate(mol.components):
                idxs.append((i,j,k))
            yield idxs


def reaction_graph_node_matchings(
    reactants: list[Pattern],
    graph_r: nx.Graph,
    graph_p: nx.Graph,
) -> list[tuple[ChemArrayNodeId, ChemArrayNodeId]]:
    """Returns a list of possible node matchings between reactants graph 
    and the products graph. Finds node matches by searching for isomorphisms
    molecule by molecule

    Args:
        reactants (list[Pattern]): input reactants
        graph_r (nx.Graph): graph corresponding to reactants
        graph_p (nx.Graph): graph corresponding to products

    Returns:
        list[tuple[ChemArrayNodeId, ChemArrayNodeId]]: list of node pairs, where the first node
        is a node in the reactants graph and the second is a node in the products graph
    """
    # We want to match nodes in graph_r to graph_p
    # However we want to do it molecule by molecule and not node by node
    # since we want to match whole molecules, between the reactants and
    # products
    node_pairs: list[tuple[ChemArrayNodeId, ChemArrayNodeId]] = []
    match_func = partial(node_matching, match_state=False, match_bond=False)

    for node_ids in reactant_idxs_iter(reactants):
        mol_subgraph = graph_r.subgraph(node_ids)
        matcher = nx.isomorphism.GraphMatcher(graph_p, mol_subgraph, match_func)
        for mappings in matcher.subgraph_isomorphisms_iter():
            for n2, n1 in mappings.items():
                node_pair = (n1, n2)
                node_pairs.append(node_pair)

    return node_pairs


def reactants_to_products_node_map(
    reactants: list[Pattern],
    graph_r: nx.Graph,
    graph_p: nx.Graph,
) -> tuple[
        dict[ChemArrayNodeId, ChemArrayNodeId],
        dict[ChemArrayNodeId, ChemArrayNodeId]
        ]:

    node_pairs = reaction_graph_node_matchings(reactants, graph_r, graph_p)

    # create a bipartite graph matching nodes in the reactants to nodes in the products
    # to find a max matching, and use the max matching to build the node map
    bipartite_graph = nx.Graph()
    top_nodes: list[tuple[int, int, int, int]] = []
    for _n1, _n2 in node_pairs:
        nb1: tuple[int, int, int, int] = (0, *_n1)
        nb2: tuple[int, int, int, int] = (1, *_n2)
        if nb1 not in bipartite_graph.nodes:
            top_nodes.append(nb1)
            data = graph_r.nodes[_n1]
            bipartite_graph.add_node(nb1, **data)
        if nb2 not in bipartite_graph.nodes:
            data = graph_p.nodes[_n2]
            bipartite_graph.add_node(nb2, **data)
        bipartite_graph.add_edge(nb1, nb2)

    max_matching = nx.algorithms.bipartite.matching.maximum_matching(
        bipartite_graph, top_nodes=top_nodes)

    # reconvert to cordinates in the reactants and products graphs
    node_map: dict[ChemArrayNodeId, ChemArrayNodeId] = {}
    node_map_reverse: dict[ChemArrayNodeId, ChemArrayNodeId] = {}
    for _n1, _n2 in max_matching.items():
        n1: ChemArrayNodeId = _n1[1:]
        n2: ChemArrayNodeId = _n2[1:]
        if _n1[0] == 0:
            node_map[n1] = n2
        else:
            node_map_reverse[n1] = n2

    return node_map, node_map_reverse


def find_state_changes(
    graph_r: nx.Graph,
    graph_p: nx.Graph,
    node_map: dict[ChemArrayNodeId, ChemArrayNodeId]
) -> dict[ChemArrayNodeId, ChemArrayNodeId]:
    """Returns a mapping between component nodes in the reactants graph and 
    component nodes in the products graph, for those components that change
    state.

    Args:
        graph_r (nx.Graph): reactants graph
        graph_p (nx.Graph): products graph
        node_map (dict[ChemArrayNodeId, ChemArrayNodeId]): node mapping

    Returns:
        dict[ChemArrayNodeId, ChemArrayNodeId]: node mapping
    """
    state_changes: dict[ChemArrayNodeId, ChemArrayNodeId] = {}
    for n1, n2 in node_map.items():
        comp1_state = graph_r.nodes[n1].get('state', None)
        comp2_state = graph_p.nodes[n2].get('state', None)

        if comp1_state != comp2_state:
            state_changes[n1] = n2
            break

    return state_changes


def find_bond_breaks(
    graph_r: nx.Graph,
    graph_p: nx.Graph,
    node_map: dict[ChemArrayNodeId, ChemArrayNodeId]
) -> tuple[
    list[ChemArrayNodeId],
    list[ChemArrayNodeId]
] | None:
    # find bonded edge on reactants graph
    for (n1, n2, bond) in graph_r.edges.data('bond', None): # type: ignore
        if not bond:
            continue

        # check if matching edge is on products graph
        n1_out = node_map.get(n1, None)
        n2_out = node_map.get(n2, None)

        # if both molecules got deleted ignore this bond
        if not n1_out and not n2_out:
            continue

        if not n1_out or not n2_out:
            raise ValueError("One of the bonded molecules was deleted.")

        # bond break
        edge_out = graph_p.edges.get((n1_out, n2_out))
        if edge_out is None:
            context_in: list[tuple[int, int, int]] = [n1, n2]
            context_out: list[tuple[int, int, int]] = [n1_out, n2_out]
            return context_in, context_out
        elif edge_out.get('bond') != bond:
            raise ValueError("Label changed from reactants to products.")


def find_bond_formations(
    graph_r: nx.Graph,
    graph_p: nx.Graph,
    node_map_reverse: dict[tuple[int, int, int], tuple[int, int, int]]
) -> tuple[
    list[ChemArrayNodeId],
    list[ChemArrayNodeId]
] | None:
    for (n1, n2, bond) in graph_p.edges.data('bond', None): # type: ignore
        if not bond:
            continue

        n1_in = node_map_reverse.get(n1, None)
        n2_in = node_map_reverse.get(n2, None)
        edge_in = graph_r.edges.get((n1_in, n2_in))
        if not n1_in or not n2_in:
            raise ValueError("One of the bonded molecules was created.")

        # bond formation
        if edge_in is None:
            context_in: list[ChemArrayNodeId] = [n1_in, n2_in]
            context_out: list[ChemArrayNodeId] = [n1, n2]
            return context_in, context_out


def find_deleted_molecules(
    reactants: list[Pattern],
    graph_r: nx.Graph,
    node_map: dict[tuple[int, int, int], tuple[int, int, int]]
):
    nodes1_g = set(node_map.keys())
    # find nodes in graph_r that are not in nodes1_g
    # these will correspond to molecules in graph_r that don't map to graph_p
    # i.e., deleted molecules
    missing_n1 = [n1 for n1 in graph_r.nodes if n1 not in nodes1_g]
    # split nodes by molecule to find which molecules are deleted
    subgraph = graph_r.subgraph(missing_n1).copy()
    for c in nx.connected_components(subgraph):
        # molecule nodes
        # crete a new pattern from each molecule and subcomponents in the connected components
        # some bonds might need to be broken
        # each connected component belongs to the same reactant
        react_idxs: set[int] = set(i for i,_,_ in c)
        assert len(react_idxs) == 1
        mol_nodes = [n1 for n1 in c if n1[2] == -1]
        molecules: list[Molecule] = [reactants[i].molecules[j].copy() for i, j, _ in mol_nodes]
        # TODO: ASSUMING THAT NO BONDS WERE BROKEN
        idx = next(iter(react_idxs))
        new_pattern = Pattern(molecules)
        if len(molecules) > 0:
            return idx, new_pattern


def find_created_molecules(
    products: list[Pattern],
    graph_p: nx.Graph,
    node_map: dict[tuple[int, int, int], tuple[int, int, int]]
):
    nodes2_g = set(node_map.values())
    # find nodes in graph_p that are not in nodes2_g
    # these will correspond to molecules in graph_p that don't map to graph_r
    # i.e., created molecules
    missing_n2 = [n2 for n2 in graph_p.nodes if n2 not in nodes2_g]

    subgraph = graph_p.subgraph(missing_n2).copy()
    for c in nx.connected_components(subgraph):
        # crete a new pattern from each molecule and subcomponents in the connected components
        # some bonds might need to be broken
        # each connected component belongs to the same reactant
        react_idxs: set[int] = set(i for i, _, _ in c)
        assert len(react_idxs) == 1
        mol_nodes = [n1 for n1 in c if n1[2] == -1]
        molecules: list[Molecule] = [products[i].molecules[j].copy()
                                     for i, j, _ in mol_nodes]
        # TODO: ASSUMING THAT NO BONDS WERE BROKEN
        new_pattern = Pattern(molecules)
        if len(molecules) > 0:
            return new_pattern


def find_bridged_vol_mol_transport(
    graph_r: nx.Graph,
    graph_p: nx.Graph,
    node_map: dict[ChemArrayNodeId, ChemArrayNodeId]
):
    # individual molecules moved from one volume to another volume through a
    # bridging surface
    # example: A(x)@Cyt -> A(x)@Nuc k

    # identify matching individual molecules
    node_groups_r: defaultdict[int,
                               set[ChemArrayNodeId]
                               ] = defaultdict(set)
    node_groups_p: defaultdict[int,
                               set[ChemArrayNodeId]
                               ] = defaultdict(set)
    for n1, n2 in node_map.items():
        node_groups_r[n1[0]].add(n1)
        node_groups_p[n2[0]].add(n2)

    for i, nodes_r in node_groups_r.items():
        n_mol = set(node[1] for node in nodes_r)


def find_transformation(
    curr_reactants: list[Pattern],
    products: list[Pattern],
    curr_graph_r: nx.Graph,
    graph_p: nx.Graph,
    node_map: dict[ChemArrayNodeId, ChemArrayNodeId],
    node_map_reverse: dict[ChemArrayNodeId, ChemArrayNodeId]
) -> list[ReactionTransformation]:

    transformations: list[ReactionTransformation] = []

    # find state changes
    state_changes = find_state_changes(
        curr_graph_r, graph_p, node_map)
    for n1, n2 in state_changes.items():
        state_change_action = ChangeStateAction(curr_reactants,
                                                curr_graph_r, graph_p, [n1], [n2])
        transformations.append(state_change_action)

    if len(transformations) != 0:
        return transformations

    # find bond breaks
    bond_breaks = find_bond_breaks(curr_graph_r, graph_p, node_map)
    if bond_breaks:
        bond_break_action = BreakBondTransform(curr_reactants, curr_graph_r, graph_p,
                                               bond_breaks[0], bond_breaks[1])
        transformations.append(bond_break_action)

    if len(transformations) != 0:
        return transformations

    # find bond formations
    bond_forms = find_bond_formations(curr_graph_r, graph_p, node_map_reverse)
    if bond_forms:
        bond_forms_action = FormBondTransform(curr_reactants, curr_graph_r, graph_p,
                                              bond_forms[0], bond_forms[1])
        transformations.append(bond_forms_action)

    # find deleted molecules
    res = find_deleted_molecules(curr_reactants, curr_graph_r, node_map)
    if res:
        idx, _ = res
        delete_molecule_action = DestroyMoleculeAction(
            curr_reactants, curr_graph_r, idx)
        transformations.append(delete_molecule_action)

    # find created molecules
    created_pattern = find_created_molecules(products, graph_p, node_map)
    if created_pattern:
        create_molecule_action = CreateMoleculeAction(
            curr_reactants, curr_graph_r, created_pattern)
        transformations.append(create_molecule_action)

    return transformations


def decompose_reaction(reaction: ReactionRule):
    curr_graph_r = reaction.reactants_graph
    graph_p = reaction.products_graph
    curr_reactants = reaction.reactants
    products = reaction.products

    transformations: list[ReactionTransformation] = []

    while not compare_chemical_array_graphs(curr_graph_r, graph_p):
        node_map, node_map_reverse = reactants_to_products_node_map(curr_reactants,
                                                                    curr_graph_r, graph_p)
        try:
            transformations2 = find_transformation(curr_reactants, products, 
                                                   curr_graph_r, graph_p,
                                                   node_map, node_map_reverse)
        except ValueError as exc:
            msg = f"Could not decompose the reaction:\n\t{reaction}"
            raise ValueError() from exc
        if len(transformations2) == 0:
            msg = (f"Could not decompose the reaction '{reaction}' " +
                   "into basic transformations.")
            raise ValueError(msg)

        transformations.extend(transformations2)
        curr_reactants = next(apply_transforms(curr_reactants, transformations2))
        curr_graph_r = build_chemical_array_graph(curr_reactants)

    return transformations
