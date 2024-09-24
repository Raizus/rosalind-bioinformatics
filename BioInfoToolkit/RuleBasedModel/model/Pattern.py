from collections import Counter, defaultdict
from itertools import product
from typing import Any, TypedDict
import networkx as nx
import graphviz
from BioInfoToolkit.RuleBasedModel.model.Component import Component, components_gen, sort_components
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Parsers import MoleculeDict, parse_pattern, parse_molecule


def components_all_equal(components: list[Component]):
    if len(components) <= 1:
        return True
    return all(components[0].matches_component(comp2) for comp2 in components[1:])


class Molecule:
    _name: str
    _components: list[Component]
    _components_counts: Counter[str]
    _compartment: str | None = None

    def __init__(self, name: str,
                 components: list[Component]) -> None:
        self._name = name
        self._components_counts = Counter(
            [component._name for component in components])
        self._compartment = None
        self._components = components

        for comp_name in self._components_counts.keys():
            comps = [comp for comp in self._components if comp._name == comp_name]
            valid = components_all_equal(comps)
            if not valid:
                msg = (f"At least one component with name {comp_name} " +
                "does not match other components with the same name.")
                raise ValueError(msg)

    @property
    def name(self):
        return self._name

    @property
    def components(self):
        return self._components

    @property
    def components_counts(self):
        return self._components_counts

    @property
    def compartment(self):
        return self._compartment

    def update_component_state(self, idx: int, new_state: str):
        if not 0 <= idx < len(self.components):
            return False

        component = self.components[idx]
        try:
            component.state = new_state
        except ValueError as _:
            return False

        return True

    @classmethod
    def from_dict(cls, parsed: MoleculeDict, molecule_types: dict[str, MoleculeType], add_remaining: bool = True):
        molecule_name = parsed["name"]
        parsed_components = parsed["components"]

        components: list[Component] = []
        all_components: list[Component] = []
        components_count: defaultdict[str, int] = defaultdict(int)

        # match to declared molecule
        if molecule_name not in molecule_types:
            raise ValueError(
                f"Molecule with name {molecule_name} was not declared.")

        # validate each component
        molecule_type = molecule_types[molecule_name]
        for parsed_component in parsed_components:
            comp_name = parsed_component['name']
            state = parsed_component['state']
            bond = parsed_component['bond']

            # check component name in molecules
            if comp_name not in molecule_type.components:
                raise ValueError(
                    f"Molecule {molecule_name} does not have a component named {comp_name}.")

            type_component = molecule_type.components[comp_name]
            states = type_component.states

            # create and validate component
            component = Component(comp_name, states, state, bond)
            if not component.matches_component(type_component):
                raise ValueError(
                    f"Component {component} is invalid, does not match component type {type_component}.")

            # check component count does not go over the maximum allowed
            max_count = molecule_type.components_counts[comp_name]
            if components_count[comp_name] >= max_count:
                raise ValueError(
                    f"Molecule {molecule_name} can only have at most {max_count} components with name {comp_name}.")

            components.append(component)
            all_components.append(component)
            components_count[comp_name] += 1

        # add remaining components
        if add_remaining:
            counts_diff: defaultdict[str, int] = defaultdict(int)

            for key in set(molecule_type.components_counts) | set(components_count):
                counts_diff[key] = molecule_type.components_counts[key] - \
                    components_count[key]

            for c_name, count in counts_diff.items():
                if count == 0:
                    continue

                states = molecule_type.components[c_name].states
                for _ in range(count):
                    comp = Component(c_name, states)
                    all_components.append(comp)
                components_count[c_name] += count

            assert components_count == molecule_type.components_counts

        molecule = Molecule(molecule_name, components)
        return molecule

    @classmethod
    def from_declaration(cls, declaration: str, molecules: dict[str, MoleculeType]):

        # parse string declaration
        parsed = parse_molecule(declaration)
        if not parsed:
            raise ValueError(f"Invalid molecule declaration: {declaration}")

        return cls.from_dict(parsed, molecules)

    def __repr__(self) -> str:
        # sort by name
        sorted_components = sort_components(self.components)
        component_str = ','.join(str(component)
                                 for component in sorted_components)
        out = f"{self._name}({component_str})"
        return out

    def is_bonded(self) -> bool:
        return any(component.is_bonded() for component in self.components)

    def match(self, other: Any, match_state: bool = False) -> bool:
        if not isinstance(other, Molecule):
            return False

        # must have the same name
        if self.name != other.name:
            return False

        # must have the same component counts
        if self.components_counts != other.components_counts:
            return False

        # must have same components
        for comp_name in self._components_counts.keys():
            comp_this = [
                comp for comp in self.components if comp.name == comp_name]
            comp_other = [
                comp for comp in self.components if comp.name == comp_name]
            comps_match = components_all_equal(comp_this) and \
                components_all_equal(comp_other) and \
                components_all_equal([comp_this[0], comp_other[0]])
            if not comps_match:
                return False

            if match_state:
                state_counts_this = Counter([comp.state for comp in comp_this])
                state_counts_other = Counter(
                    [comp.state for comp in comp_other])
                if state_counts_this != state_counts_other:
                    return False

        return True

    def copy(self) -> "Molecule":
        new_components = [component.copy() for component in self.components]
        new = Molecule(self._name, new_components)
        return new

    def is_specie(self) -> bool:
        return all(component.is_stateless() or (component.state in component.state)
                   for component in self.components)

    def is_pattern(self) -> bool:
        return not self.is_specie()

    def generate_species(self):
        if self.is_specie():
            yield self
            return

        components_iters = [components_gen(component)
                            for component in self.components]
        for comps in product(*components_iters):
            specie = Molecule(self._name, list(comps))
            yield specie


def generate_species(molecule: MoleculeType):
    name = molecule.name

    molecule_type_components = [component for name, component in molecule.components.items(
    ) for _ in range(molecule.components_counts[name])]
    components_iters = [
        components_gen(component) for component in molecule_type_components]

    for aux2 in product(*components_iters):
        reagent = Molecule(name, list(aux2))
        yield reagent


#     def generate_species(self):
#         if self.is_specie():
#             yield self
#             return

#         molecules_iters = [molecule.generate_species() for molecule in self.parts]
#         for molecules in product(*molecules_iters):
#             specie = ComplexReactant(list(molecules))
#             yield specie

#             # TODO: Account for multiplicity
#             # find bonded components with multiplicity greater than 1
#             # generate all possible combinations

#     def number_molecules(self) -> int:
#         return len(self.parts)


# class NodeLabelDict(TypedDict):
#     molecule_name: str
#     component_name: str | None
#     state: str | None
#     bond: str | None


class Pattern:
    molecules: list[Molecule]
    _bonds: defaultdict[str, list[tuple[int, int]]]
    _graph: nx.Graph

    def __init__(self, molecules: list[Molecule]) -> None:
        self.molecules = molecules
        self._bonds = defaultdict(list)

        graph = nx.Graph()
        self._graph = graph
        # create graph
        for i, molecule in enumerate(self.molecules):
            parent_id = (i, -1)
            label = (molecule._name, None, None, None)
            graph.add_node(parent_id, label=label, node_id=parent_id)
            for j, component in enumerate(molecule.components):
                child_id = (i, j)
                label = (molecule._name, component.name,
                            component.state, component.bond)
                graph.add_node(child_id, label=label, node_id=child_id)
                graph.add_edge(parent_id, child_id)
                bond = component.bond
                if len(bond) and bond != '?':
                    self._bonds[bond].append(child_id)

        for bond_id, nodes in self._bonds.items():
            if len(nodes) == 2 and bond_id != '+':
                node1, node2 = nodes
                graph.add_edge(node1, node2, bond=bond_id)
            else:
                msg = f"There are bonds ({bond_id}) with a number of nodes different than 2 ({len(nodes)})."
                print(msg)

    @classmethod
    def from_dict(cls,
                  parsed: list[MoleculeDict],
                  molecule_types: dict[str, MoleculeType]):
        parts = []
        for parsed_reagent in parsed:
            part = Molecule.from_dict(parsed_reagent, molecule_types)
            parts.append(part)
        return Pattern(parts)

    @classmethod
    def from_declaration(cls,
                         declaration: str,
                         molecule_types: dict[str, MoleculeType]):
        parsed = parse_pattern(declaration)
        reactant = cls.from_dict(parsed, molecule_types)
        return reactant

    @property
    def graph(self):
        return self._graph

    def __eq__(self, other: object) -> bool:
        """Compares to patterns. If two different pattern graphs are isomorphic, then
        the patterns are considered equal. Nodes match if they are the same molecule, same_component, same_component state and bonding. The bonding label may be different but the edges have to match (the edge label may not) (obviously, or it wouldn't be isomorphic)

        Args:
            other (object): _description_

        Returns:
            bool: _description_
        """
        if not isinstance(other, Pattern):
            return False

        counts1 = self.molecule_counts()
        counts2 = other.molecule_counts()
        # check for the molecules
        if counts1 != counts2:
            return False

        def node_matching(n1: Any, n2: Any) -> bool:
            color_id1: tuple[str, str | None, str | None, str | None] = n1['label']
            color_id2: tuple[str, str | None, str | None, str | None] = n2['label']
            mol1_name, comp1_name, comp1_state, comp1_bond = color_id1
            mol2_name, comp2_name, comp2_state, comp2_bond = color_id2

            component_match = mol1_name == mol2_name and comp1_name == comp2_name
            state_match = comp1_state == comp2_state

            # check bonds
            comp1_is_bonded = comp1_bond is not None and len(comp1_bond) > 0 and comp1_bond != '?'
            comp2_is_bonded = comp2_bond is not None and len(comp2_bond) > 0 and comp2_bond != '?'

            bond_match = comp1_is_bonded == comp2_is_bonded
            if comp1_bond == '?' or comp2_bond == '?':
                bond_match = True
            elif comp1_bond == '+':
                bond_match = comp2_bond is not None and len(comp2_bond) > 0


            full_match = component_match and state_match and bond_match
            return full_match

        is_iso = nx.isomorphism.is_isomorphic(self.graph, other.graph, node_matching)
        return is_iso

    def molecule_counts(self) -> Counter[str]:
        counts = Counter(molecule.name for molecule in self.molecules)
        return counts

    def update_component_state(self, pos: tuple[int, int], new_state: str):
        i, j = pos
        if not 0 <= i < len(self.molecules):
            return False

        molecule = self.molecules[i]
        if not 0 <= j < len(molecule.components):
            return False

        component = molecule.components[j]
        try:
            component.state = new_state
            # need to update the graph
            label: tuple[str, str | None, str | None, str |
                         None] = self.graph.nodes[(i, j)]['label']
            new_label = (label[0], label[1], new_state, label[3])
            nx.set_node_attributes(self.graph, {(i, j): new_label}, 'label')
        except ValueError as _:
            return False

        return True

    def number_of_molecules(self) -> int:
        return len(self.molecules)

    def is_connected(self) -> bool:
        """Returns True if the graph of the pattern is connected, False otherwise

        Returns:
            bool: connected
        """
        connected = nx.is_connected(self._graph)
        return connected

    def __repr__(self) -> str:
        out = '.'.join(str(part) for part in self.molecules)
        return out

    def copy(self) -> "Pattern":
        molecules = self.molecules
        molecules = [m.copy() for m in molecules]
        new = Pattern(molecules)
        return new

    def draw_graph(self, filename='pattern_graph'):
        dot = graphviz.Graph(comment='Molecule Pattern Graph', format='png')

        # Add nodes for molecules and components
        for node_id, node_data in self._graph.nodes(data=True):
            label_data = node_data['label']
            molecule_name = label_data[0]
            component_name = label_data[1]
            state = label_data[2]
            bond = label_data[3]
            label = ""
            shape = ""
            style = ""

            # Molecule node
            if component_name is None:
                label = f"({molecule_name}, {node_id[0]})"
                shape = 'rectangle'
            else:
                # Component node
                molecule_idx, component_idx = node_id
                label = f"({component_name}, {molecule_idx}, {component_idx})"
                if state:
                    label += f"~{state}"
                shape = 'ellipse'

                # Determine border style
                if bond == '?':
                    style = 'dashed'
                elif bond:
                    style = 'dotted'
                else:
                    style = 'solid'

            dot.node(str(node_id), label=label, shape=shape, style=style)

        # Add edges between molecules and components
        for edge_start, edge_end, edge_data in self._graph.edges(data=True):
            bond = edge_data.get('bond', None)
            style = 'solid'
            label = ""

            # If it's a bond between two components
            if bond:
                style = 'dashed'
                label = bond

            dot.edge(str(edge_start), str(edge_end), style=style, label=label)

        # Render graph to file
        dot.render(filename)


def node_matching_func(n1: Any, n2: Any) -> bool:
    color_id1: tuple[str, str | None, str | None, str | None] = n1['label']
    color_id2: tuple[str, str | None, str | None, str | None] = n2['label']
    mol1_name, comp1_name, comp1_state, comp1_bond = color_id1
    mol2_name, comp2_name, comp2_state, comp2_bond = color_id2

    component_match = mol1_name == mol2_name and comp1_name == comp2_name
    state_match = comp1_state == comp2_state
    if comp2_state is not None and len(comp2_state) == 0:
        state_match = True

    # check bonds
    bond_match = comp1_bond == comp2_bond
    if comp2_bond == '?':
        bond_match = True
    elif comp2_bond == '+':
        bond_match = comp1_bond is not None and len(comp1_bond) > 0

    full_match = component_match and state_match and bond_match
    return full_match


def match_patterns(pattern1: Pattern, pattern2: Pattern, count: bool = False) -> int:
    """Checks if the pattern1 graph is isomorphic to a subgraph of the pattern2 graph.
    For species-observables, count should be False and for molecules-observables it should be True.
    If count is False returns 1 if there's a match and 0 otherwise.
    If count is True returns a count of all subgraph isomorphisms with distinct nodes.
    """
    counts1 = pattern1.molecule_counts()
    counts2 = pattern2.molecule_counts()

    # check for the molecules
    for name, count1 in counts1.items():
        count2 = counts2.get(name, 0)
        if count2 < count1:
            print(f"{pattern2} needs to have at least {count1} '{name}' molecules")
            return 0

    # check for subgraph isomorphism
    matcher = nx.isomorphism.GraphMatcher(
        pattern2.graph, pattern1.graph, node_matching_func)
    num_matches = 0
    if count:
        nodes_set: list[set[tuple[int, int]]] = []
        num_matches = 0
        for g in matcher.subgraph_isomorphisms_iter():
            nodes = set(n for n in g.keys())
            if nodes not in nodes_set:
                nodes_set.append(nodes)
                num_matches += 1
        # num_matches = len([1 for _ in matcher.subgraph_isomorphisms_iter()])
    else:
        num_matches = int(matcher.subgraph_is_isomorphic())
    return num_matches


# if __name__ == "__main__":
#     pattern_str = "A(x,y!0).B(p!0,q~q1!1).C(r~r2!1,s!2).A(x!2,y)"
#     molecule_types = {
#         "A": MoleculeType("A(x,y)"),
#         "B": MoleculeType("B(p,q~q1~q2)"),
#         "C": MoleculeType("C(r~r1~r2,s)"),
#     }

#     pattern = Pattern.from_declaration(pattern_str, molecule_types)
#     pattern.draw_graph()
