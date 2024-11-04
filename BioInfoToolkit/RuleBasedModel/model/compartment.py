"""
Module for Compartment and Compartments classes, to build compartmentalized models.
"""

from typing import Collection, OrderedDict
import networkx as nx

from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_compartment
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import CompartmentDict


class Compartment:
    """Class representing a cell compartment"""
    name: str
    dimensions: int
    volume: str
    enclosing_compartment: str | None
    comment: str

    def __init__(
        self,
        compartment_dict: CompartmentDict
    ) -> None:
        self.name = compartment_dict["name"]
        self.dimensions = compartment_dict["dimensions"]
        self.volume = compartment_dict["volume"]
        self.enclosing_compartment = compartment_dict["enclosing_compartment"]
        self.comment = compartment_dict["comment"]

    @classmethod
    def from_declaration(cls, declaration: str) -> "Compartment":
        """Creates a Compartment from the string declaring the compartment

        Args:
            declaration (str): the string declaration of the compartment

        Returns:
            compartment (Compartment): compartment
        """
        parsed = parse_compartment(declaration)
        compartment = Compartment(parsed)
        return compartment

    def is_enclosed(self) -> bool:
        """Returns True if the compartment is enclosed (has a parent),
        False otherwise

        Returns:
            enclosed (bool): bool indicating if compartment is enclosed
        """
        return bool(self.enclosing_compartment)

    def __repr__(self) -> str:
        out = f"{self.name} {self.dimensions} {self.volume} {self.enclosing_compartment}"
        return out


def find_path(
    compartment_tree: nx.DiGraph,
    nodes: Collection[str]
) -> list[str] | None:
    """Finds a path in the compartment tree that 
    only contains all of the given nodes

    Args:
        compartment_tree (nx.DiGraph): compartment tree
        nodes (Collection[str]): collection of nodes that form the path

    Returns:
        path (list[str] | None): list of nodes forming the path, if it exists
    """

    def path_to_target(root: str, target: str):
        # Helper function to find the path from a root to a target node
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
    """Class to represent the compartmentalization structure of the model.
    """
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
        """Returns True if the two compartments are adjacent, False otherwise.

        Args:
            comp1 (str): first compartment
            comp2 (str): second compartment

        Returns:
            adjacent (bool): bool indicating if the given compartments are adjacent
        """
        tree = self.compartment_tree
        if tree.has_edge(comp1, comp2) or tree.has_edge(comp2, comp1):
            return True
        return False

    def get_aggregate(self, comps: Collection[str]) -> str | None:
        """Given a collection of compartments, returns the resulting aggregate compartment 
        if it exists, and returns None if the collection is empty.
        Raises a value error if the compartment collection would result in a topologically
        inconsistent species / pattern

        Args:
            comps (Collection[str]): collection of compartments

        Raises:
            ValueError: error raised if compartment collection would result in a 
            topologically inconsistent species / pattern

        Returns:
            aggregate_compartment (str | None): the aggregate compartment
        """
        # check all compartments are valid
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
            msg = "The compartments spaned by a specie must be adjacent to each other."
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
        """Returns True if there is at least one compartment, False otherwise.

        Returns:
            compartmentalized (bool): bool indicating if there are any compartments defined
        """
        return len(self.items) > 0

    def has_compartment(self, compartment: str) -> bool:
        """Returns True if the the given compartment exists, False otherwise

        Args:
            compartment (str): compartment name

        Returns:
            compartment_exists (bool): bool indicating if compartment exists
        """
        if self.items.get(compartment):
            return True
        return False


def is_valid_surface(
    compartment: str,
    comps_dict: dict[str, Compartment],
    comp_tree: nx.DiGraph
) -> bool:
    """Checks if a compartment is a valid surface

    Args:
        compartment (str): compartment name
        comps_dict (dict[str, Compartment]): compartments dictionary
        comp_tree (nx.DiGraph): compartments tree

    Returns:
        valid_surface (bool): __description__
    """
    dim = comps_dict[compartment].dimensions

    if dim == 2:
        # is a surface, must only have 1 volume child
        count = 0
        for child in comp_tree.successors(compartment):
            count += 1
            child_dim = comps_dict[child].dimensions
            if child_dim != 3:
                return False
        if count != 1:
            return False

        # is a surface, must only have 1 volume parent
        count = 0
        for parent in comp_tree.predecessors(compartment):
            count += 1
            parent_dim = comps_dict[parent].dimensions
            if parent_dim != 3:
                return False
        if count != 1:
            return False

        return True

    return False


def is_valid_volume(
    compartment: str,
    comps_dict: dict[str, Compartment],
    comp_tree: nx.DiGraph
) -> bool:
    """Checks if a compartment is a valid volume

    Args:
        compartment (str): compartment name
        comps_dict (dict[str, Compartment]): compartments dictionary
        comp_tree (nx.DiGraph): compartments tree

    Returns:
        valid_volume (bool): __description__
    """
    dim = comps_dict[compartment].dimensions

    if dim == 3:
        # is a volume, can have any number of children,
        # but they must all be surfaces
        for child in comp_tree.successors(compartment):
            child_dim = comps_dict[child].dimensions
            if child_dim != 2:
                return False

        # at most 1 parent that must be a surface
        count = 0
        for parent in comp_tree.predecessors(compartment):
            count += 1
            parent_dim = comps_dict[parent].dimensions
            if parent_dim != 2:
                return False
        if count > 1:
            return False

        return True

    return False
