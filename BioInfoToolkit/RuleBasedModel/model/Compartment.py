
from typing import Collection, OrderedDict
import networkx as nx

from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_compartment
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import CompartmentDict


class Compartment:
    name: str
    dimensions: int
    volume: str
    enclosing_compartment: str | None
    comment: str

    def __init__(self, name: str, dimensions: int,
                 volume: str, enclosing_compartment: str | None,
                 comment: str) -> None:
        self.name = name
        self.dimensions = dimensions
        self.volume = volume
        self.enclosing_compartment = enclosing_compartment
        self.comment = comment

    @classmethod
    def from_dict(cls, parsed: CompartmentDict) -> "Compartment":
        name = parsed["name"]
        dimensions = parsed["dimensions"]
        volume = parsed["volume"]
        enclosing_compartment = parsed["enclosing_compartment"]
        comment = parsed["comment"]

        compartment = Compartment(name, dimensions, volume, enclosing_compartment, comment)
        return compartment

    @classmethod
    def from_declaration(cls, declaration: str) -> "Compartment":
        parsed = parse_compartment(declaration)
        compartment = cls.from_dict(parsed)
        return compartment

    def is_enclosed(self) -> bool:
        return bool(self.enclosing_compartment)

    def __repr__(self) -> str:
        out = f"{self.name} {self.dimensions} {self.volume} {self.enclosing_compartment}"
        return out


def find_path(
    compartment_tree: nx.DiGraph,
    nodes: Collection[str]
) -> list[str] | None:
    # Helper function to find the path from a root to a target node
    
    def path_to_target(root, target):
        try:
            return nx.shortest_path(compartment_tree, source=root, target=target)
        except nx.NetworkXNoPath:
            return None

    if compartment_tree.number_of_nodes() == 0:
        return None

    # Find the lowest common ancestor (LCA) of the nodes in the tree
    root = next(
        n for n in compartment_tree if compartment_tree.in_degree(n) == 0)
    common_ancestor = nx.lowest_common_ancestor(compartment_tree, root, *nodes)

    if common_ancestor is None:
        return None

    # Find paths from the LCA to each node in the collection
    paths = []
    for node in nodes:
        path = path_to_target(common_ancestor, node)
        if path is None:
            return None
        paths.append(path)

    # Construct a single path that passes through all nodes in order
    final_path = []
    visited = set()
    for path in paths:
        for node in path:
            if node not in visited:
                final_path.append(node)
                visited.add(node)

    # Verify final path contains only the specified nodes
    if set(final_path) == set(nodes):
        return final_path

    return None


class Compartments:
    items: OrderedDict[str, Compartment]
    compartment_tree: nx.DiGraph
    root: str | None

    def __init__(self,
        items: OrderedDict[str, Compartment],
        compartment_tree: nx.DiGraph,
        root: str | None
    ) -> None:
        self.items = items
        self.root = root
        self.compartment_tree = compartment_tree

    def are_adjacent(self, comp1: str, comp2: str) -> bool:
        tree = self.compartment_tree
        if tree.has_edge(comp1, comp2) or tree.has_edge(comp2, comp1):
            return True
        return False

    def get_aggregate(self, comps: Collection[str]):
        # check all components are valid
        invalid_comps = set(comp for comp in comps if comp not in self.items)
        if len(invalid_comps) > 0:
            msg = f"Invalid compartment(s): {invalid_comps}."
            raise ValueError(msg)

        if len(comps) == 0:
            return None

        if len(comps) == 1:
            aggregate = next(iter(comps))
            return aggregate

        if len(comps) == 2:
            comp1, comp2 = list(comps)
            if self.are_adjacent(comp1, comp2):
                aggregate = comp1 if self.items[comp1].dimensions == 2 else comp2
                return aggregate
            msg = "The components spaned by a specie must be adjacent to each other."
            raise ValueError(msg)

        if len(comps) == 3:
            path = find_path(self.compartment_tree, comps)
            if path is None:
                msg = ("A specie can only span at most 3 compartments, "
                       "1 surface and 2 adjacent volumes.")
                raise ValueError(msg)

            surfaces: list[str] = []
            volumes: list[str] = []
            for comp in comps:
                if self.items[comp].dimensions == 2:
                    surfaces.append(comp)
                else:
                    volumes.append(comp)

            if len(surfaces) == 1:
                aggregate = surfaces[0]
                return aggregate

            msg = "A specie can only span at most 1 surface."
            raise ValueError(msg)

        msg = "A specie can only span at most 3 compartments, 1 surface and 2 adjacent volumes."
        raise ValueError(msg)

    def is_compartmentalized(self) -> bool:
        return len(self.items) > 0

    def has_compartment(self, compartment: str) -> bool:
        if self.items.get(compartment):
            return True
        return False

def is_valid_surface(
    compartment: str,
    comps_dict: dict[str, Compartment],
    comp_graph: nx.DiGraph
):
    dim = comps_dict[compartment].dimensions

    if dim == 2:
        # is a surface, must only have 1 volume child
        count = 0
        for child in comp_graph.successors(compartment):
            count += 1
            child_dim = comps_dict[child].dimensions
            if child_dim != 3:
                return False
        if count != 1:
            return False
        return True
    return False


def is_valid_volume(
    compartment: str,
    comps_dict: dict[str, Compartment],
    comp_graph: nx.DiGraph
    ):
    dim = comps_dict[compartment].dimensions

    if dim == 3:
        # is a volume, can have any numbe of children,
        # but they must all be surfaces
        for child in comp_graph.successors(compartment):
            child_dim = comps_dict[child].dimensions
            if child_dim != 2:
                return False
        return True
    return False
