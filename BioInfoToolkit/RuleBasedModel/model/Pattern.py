from collections import Counter, defaultdict
from itertools import product
from typing import Any
import networkx as nx
import graphviz
from BioInfoToolkit.RuleBasedModel.model.Compartment import Compartments
from BioInfoToolkit.RuleBasedModel.model.Component import Component, components_all_equal, components_gen, sort_components
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_pattern, parse_molecule
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import MoleculeDict, PatternDict
from BioInfoToolkit.RuleBasedModel.utils.utls import apply_inequality


class Molecule:
    _name: str
    _components: list[Component]
    _components_counts: Counter[str]
    _compartment: str | None = None

    def __init__(self, name: str,
                 components: list[Component],
                 compartment: str | None = None) -> None:
        self._name = name
        sorted_comps = sort_components(components)
        self._components = sorted_comps
        self._components_counts = Counter(
            [component._name for component in sorted_comps])
        self._compartment = compartment

        for comp_name in self._components_counts.keys():
            comps = [comp for comp in self._components if comp._name == comp_name]
            valid = components_all_equal(comps)
            if not valid:
                msg = (f"At least one component with name {comp_name} " +
                "does not match other components with the same name.")
                raise ValueError(msg)

    def validate(self, molecule_types: dict[str, MoleculeType] | None) -> bool:
        name = self.name
        if molecule_types and name not in molecule_types:
            return False

        molecule_type = molecule_types[name] if molecule_types else None
        if molecule_type:
            for comp in self.components:
                if comp.name not in molecule_type.components:
                    return False

                type_component = molecule_type.components[comp.name]
                if not comp.matches_component(type_component):
                    return False

            molecule_comp_counts = self.components_counts
            type_comp_counts = molecule_type.components_counts

            for comp_name, count in molecule_comp_counts.items():
                if count > type_comp_counts.get(comp_name, 0):
                    return False

        return True

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

    @compartment.setter
    def compartment(self, compartment: str):
        self._compartment = compartment

    @classmethod
    def from_dict(cls, parsed: MoleculeDict,
                  molecule_types: dict[str, MoleculeType] | None):
        molecule_name = parsed["name"]
        parsed_components = parsed["components"]
        compartment = parsed['compartment']

        components: list[Component] = []
        all_components: list[Component] = []
        components_count: defaultdict[str, int] = defaultdict(int)

        # match to declared molecule
        if molecule_types and molecule_name not in molecule_types:
            raise ValueError(
                f"Molecule with name {molecule_name} was not declared.")

        # validate each component
        molecule_type = molecule_types[molecule_name] if molecule_types else None
        for parsed_component in parsed_components:
            comp_name = parsed_component['name']
            state = parsed_component['state']
            bond = parsed_component['bond']

            # check component name in molecules
            if molecule_type and comp_name not in molecule_type.components:
                msg = f"Molecule {molecule_name} does not have a component named {comp_name}."
                raise ValueError(msg)

            # create component
            if molecule_type:
                type_component = molecule_type.components[comp_name]
                states: set[str] | None = type_component.states

                # create and validate component
                component = Component(comp_name, states, state, bond)
                if not component.matches_component(type_component):
                    msg = (f"Component {component} is invalid, " +
                           f"does not match component type {type_component}.")
                    raise ValueError(msg)
            else:
                states: set[str] | None = None
                component = Component(comp_name, states, state, bond)

            # check component count does not go over the maximum allowed
            if molecule_type:
                max_count = molecule_type.components_counts[comp_name]
                if components_count[comp_name] >= max_count:
                    msg = (f"Molecule {molecule_name} can only have at most {max_count} " +
                        f"components with name {comp_name}.")
                    raise ValueError(msg)

            components.append(component)
            all_components.append(component)
            components_count[comp_name] += 1

        molecule = Molecule(molecule_name, components, compartment)
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
        if self.compartment:
            out += f"@{self.compartment}"
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

    def change_state(self, idx: int, state: str) -> bool:
        if 0 <= idx < len(self.components):
            comp = self.components[idx]
            comp.state = state
            return True
        return False


def sort_molecules(molecules: list[Molecule]) -> list[Molecule]:
    def sort_fun(x: Molecule):
        comps_str = ','.join(str(comp) for comp in x.components)
        return (x.name, comps_str)

    sorted_mol = sorted(molecules, key=sort_fun)
    return sorted_mol

def generate_species(molecule: MoleculeType):
    name = molecule.name

    molecule_type_components = [component for name, component in molecule.components.items(
    ) for _ in range(molecule.components_counts[name])]
    components_iters = [
        components_gen(component) for component in molecule_type_components]

    for aux2 in product(*components_iters):
        reagent = Molecule(name, list(aux2))
        yield reagent


class Pattern:
    molecules: list[Molecule]
    aggregate_compartment: str | None
    _bonds: defaultdict[str, list[tuple[int, int]]]
    _graph: nx.Graph

    def __init__(self, molecules: list[Molecule],
                 aggregate_compartment: str | None = None) -> None:
        self.molecules = sort_molecules(molecules)
        self._bonds = defaultdict(list)
        self.aggregate_compartment = aggregate_compartment

        graph = nx.Graph()
        self._graph = graph
        # create graph
        for i, molecule in enumerate(self.molecules):
            parent_id = (i, -1)
            graph.add_node(parent_id,
                           molecule_name=molecule.name,
                           node_id=parent_id,
                           compartment=molecule.compartment)
            
            for j, component in enumerate(molecule.components):
                child_id = (i, j)
                graph.add_node(child_id, molecule_name=molecule.name,
                               comp_name=component.name, states=component.states,
                               state=component.state, bond=component.bond,
                               node_id=child_id)
                graph.add_edge(parent_id, child_id)
                bond = component.bond
                if len(bond) and bond not in ('?', '+'):
                    self._bonds[bond].append(child_id)

        for bond_id, nodes in self._bonds.items():
            if len(nodes) == 2:
                node1, node2 = nodes
                graph.add_edge(node1, node2, bond=bond_id)
            else:
                msg = (f"There are bonds ({bond_id}) with a number " +
                "of nodes different than 2 ({len(nodes)}).")
                print(msg)

    @classmethod
    def from_dict(cls,
                  parsed: PatternDict,
                  molecule_types: dict[str, MoleculeType] | None):
        parts = []
        for parsed_mol in parsed['molecules']:
            molecule = Molecule.from_dict(parsed_mol, molecule_types)
            parts.append(molecule)
        aggregate_compartment = parsed['aggregate_compartment']
        return Pattern(parts, aggregate_compartment)

    @classmethod
    def from_declaration(cls,
                         declaration: str,
                         molecule_types: dict[str, MoleculeType]):
        parsed = parse_pattern(declaration)
        reactant = cls.from_dict(parsed, molecule_types)
        return reactant

    def validate(
        self,
        molecule_types: dict[str, MoleculeType] | None,
        compartments: Compartments | None = None
    ) -> bool:
        for mol in self.molecules:
            if not mol.validate(molecule_types):
                return False

        if compartments and not self.is_topologically_consistent(compartments):
            return False

        return True

    @property
    def graph(self):
        return self._graph

    @property
    def bonds(self):
        return self._bonds

    def __eq__(self, other: object) -> bool:
        """Compares to patterns. If two different pattern graphs are isomorphic, then
        the patterns are considered equal. Nodes match if they are the same molecule,
        same_component, same_component state and bonding. The bonding label may be
        different but the edges have to match (the edge label may not).

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
            mol1_name = n1.get('molecule_name')
            comp1_name = n1.get('comp_name', None)
            comp1_state = n1.get('state', None)
            comp1_bond = n1.get('bond', None)

            mol2_name = n2.get('molecule_name')
            comp2_name = n2.get('comp_name', None)
            comp2_state = n2.get('state', None)
            comp2_bond = n2.get('bond', None)

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
        out = '.'.join(str(mol) for mol in self.molecules)
        if self.aggregate_compartment:
            out = f"@{self.aggregate_compartment}:{out}"
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
            molecule_name: str | None = node_data.get('molecule_name')
            component_name: str | None = node_data.get('comp_name', None)
            state: str | None = node_data.get('state', None)
            bond: str | None = node_data.get('bond', None)

            label = ""
            shape = ""
            style = ""

            # Molecule node
            if component_name is None:
                label = f"({molecule_name})"
                shape = 'rectangle'
            else:
                # Component node
                label = f"({component_name}"
                if state:
                    label += f"~{state}"
                label += ')'
                shape = 'ellipse'

                # Determine border style
                if bond == '?':
                    style = 'dashed'
                elif bond:
                    style = 'dotted'
                else:
                    style = 'solid'

            label = f"{node_id}\n"+ label
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

    def break_bond(self, bond_id: str):
        n1, n2 = self.bonds[bond_id]
        g_copy = self.graph.copy()
        g_copy.remove_edge(n1, n2)

        molecules = [mol.copy() for mol in self.molecules]
        for n in (n1, n2):
            molecules[n[0]].components[n[1]].bond = ''

        # build 2 new patterns
        patterns: list[Pattern] = []
        for nodes in nx.connected_components(g_copy):
            molecule_idxs = [i for (i, j) in nodes if j == -1]
            molecules_a = [molecules[i] for i in molecule_idxs]
            new_patt = Pattern(molecules_a)
            patterns.append(new_patt)

        pattern1, pattern2 = patterns[0], patterns[1]
        return pattern1, pattern2

    def is_topologically_consistent(self, compartments: Compartments) -> bool:
        # a pattern is topologically consistent if it does not span more than one surface
        # and if it does not have bonds that need to pass through compartments
        spaned_compartments: set[str] = set()
        if self.aggregate_compartment:
            spaned_compartments.add(self.aggregate_compartment)

        for mol in self.molecules:
            if mol.compartment:
                spaned_compartments.add(mol.compartment)

        try:
            aggregate = compartments.get_aggregate(spaned_compartments)
            if self.aggregate_compartment and aggregate != self.aggregate_compartment:
                return False
        except ValueError:
            return False

        # check that bonds don't cross compartments
        # bonds must be either in the same compartment or between adjacent compartments
        for _, bond_nodes in self.bonds.items():
            if len(bond_nodes) != 2:
                return False
            mol1_idx, mol2_idx = bond_nodes[0][0], bond_nodes[1][0]
            compart1 = self.molecules[mol1_idx].compartment
            compart1 = compart1 if compart1 else self.aggregate_compartment
            compart2 = self.molecules[mol2_idx].compartment
            compart2 = compart2 if compart2 else self.aggregate_compartment
            if ((compart1 and compart2)
                and ((compart1 != compart2)
                     and not compartments.are_adjacent(compart1, compart2))):
                return False
        return True


def split_pattern_into_connected_components(patt: Pattern):
    molecules = [mol.copy() for mol in patt.molecules]
    patterns: list[Pattern] = []
    for nodes in nx.connected_components(patt.graph):
        molecule_idxs = [i for (i, j) in nodes if j == -1]
        molecules_a = [molecules[i] for i in molecule_idxs]
        new_patt = Pattern(molecules_a)
        patterns.append(new_patt)
    return patterns


def form_bond(pattern1: Pattern,
              pattern2: Pattern,
              n1: tuple[int, int],
              n2: tuple[int, int]) -> Pattern:

    molecules1 = [mol.copy() for mol in pattern1.molecules]
    molecules2 = [mol.copy() for mol in pattern2.molecules]
    # we must rename the bonds to avoid collision
    bond_id = 1
    for molecules in (molecules1, molecules2):
        bond_remap: dict[str, str] = {}
        for mol in molecules:
            for comp in mol.components:
                if not (comp.bond and comp.bond not in ('+', '?')):
                    continue

                if comp.bond in bond_remap:
                    comp.bond = bond_remap[comp.bond]
                else:
                    bond_remap[comp.bond] = str(bond_id)
                    comp.bond = bond_remap[comp.bond]
                    bond_id += 1

    molecules1[n1[0]].components[n1[1]].bond = str(bond_id)
    molecules2[n2[0]].components[n2[1]].bond = str(bond_id)

    new_pattern = Pattern(molecules1+molecules2)
    return new_pattern


def node_pattern_matching_func(n1: Any, n2: Any) -> bool:
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
    if comp2_state is not None and len(comp2_state) == 0:
        state_match = True

    # check bonds
    is_bonded1 = bool(comp1_bond)
    is_bonded2 = bool(comp2_bond)
    bond_match = is_bonded1 == is_bonded2
    if comp2_bond == '?' or comp1_bond == '?':
        bond_match = True
    elif comp2_bond == '+':
        bond_match = comp1_bond is not None and len(comp1_bond) > 0
    elif comp1_bond == '+':
        bond_match = comp2_bond is not None and len(comp2_bond) > 0

    full_match = component_match and state_match and bond_match
    return full_match


def match_pattern_specie(
    pattern: Pattern,
    specie: Pattern,
    count_unique_matches: bool = False
) -> int:
    """Checks if the pattern graph is isomorphic to a subgraph of the specie graph.
    To check if a species pattern matches a general pattern, specie should be the species
    pattern and pattern1 the general pattern.
    For species-observables, count should be False and for molecules-observables it should
    be True.

    Args:
        pattern1 (Pattern): _description_
        specie (Pattern): _description_
        count_unique_matches (bool, optional): If count is False returns 1 if there's a match 
        and 0 otherwise. If count is True returns a count of all subgraph isomorphisms with 
        distinct nodes. Defaults to False. 

    Returns:
        int: _description_
    """
    counts1 = pattern.molecule_counts()
    counts2 = specie.molecule_counts()

    # check for the molecules
    for name, count1 in counts1.items():
        count2 = counts2.get(name, 0)
        if count2 < count1:
            # print(f"{specie} needs to have at least {count1} '{name}' molecules")
            return 0

    # check for subgraph isomorphism
    match_func = node_pattern_matching_func
    matcher = nx.isomorphism.GraphMatcher(
        specie.graph, pattern.graph, match_func)

    num_matches = 0
    if count_unique_matches:
        nodes_set: list[set[
            tuple[tuple[int, int], tuple[int, int]]
            ]] = []
        num_matches = 0
        for mapping in matcher.subgraph_isomorphisms_iter():
            nodes = set((n1, n2) for n1, n2 in mapping.items())
            if nodes not in nodes_set:
                nodes_set.append(nodes)
                num_matches += 1

    else:
        num_matches = int(matcher.subgraph_is_isomorphic())

    return num_matches
