import abc
from collections import defaultdict
from enum import Enum
from typing import Any, Generator
import networkx as nx

from BioInfoToolkit.RuleBasedModel.model.Pattern import Molecule, Pattern, form_bond, \
    sp_pattern_node_matching, split_pattern_into_connected_components
from BioInfoToolkit.RuleBasedModel.model.RuleModifiers import RuleModifiers

def idxs_from_reactants(reactants: list[Pattern]):
    idxs: list[tuple[int, int, int]] = []
    for i, reactant in enumerate(reactants):
        for j, molecule in enumerate(reactant.molecules):
            for k, _ in enumerate(molecule.components):
                idxs.append((i, j, k))

    return idxs


def objs_from_idx(reactants: list[Pattern], idx: tuple[int, int, int]):
    i, j, k = idx
    reactant = reactants[i]
    molecule = reactant.molecules[j]
    component = molecule.components[k]
    return reactant, molecule, component


class TransformationType(Enum):
    BREAK_BOND = "Break Bond"
    FORM_BOND = "Form Bond"
    CHANGE_COMPONENT_STATE = "Change Component State"
    CREATE_MOLECULE = "Create Molecules"
    DESTROY_MOLECULE = "Destroy Molecule"
    ADJ_COMPARTMENT_MOL_TRANSPORT = "Adjacent Compartment Molecular Transport"
    BRIDGED_VOLUME_MOL_TRANSPORT = "Bridged-Volume Molecular Transport"
    BRIDGED_SURFACE_CONNECTED_TRANSPORT = "Bridged-Surface Connected Transport"


class ReactionTransformation(abc.ABC):
    transformation: TransformationType
    # tuple with (reactant_idx, molecule_idx, component_idx)
    # list of components that participate in a tranformation
    #   - Form bond: involves 2 components in two different molecules and different reactants
    #   - Break bond: involves 2 components in two different molecules in the same reactants
    #   - Change Component State: Involves a single component

    # the reaction center here refers to the reaction rules
    # when applying the action to a given chemical array,
    # the corresponding node mapping must be provided
    reaction_center_in: list[tuple[int, int, int]]
    reaction_center_out: list[tuple[int, int, int]]
    patterns_graph: nx.Graph

    def __init__(self, graph: nx.Graph) -> None:
        self.reaction_center_in = []
        self.reaction_center_out = []
        self.patterns_graph = graph

    @abc.abstractmethod
    def apply(
        self,
        reactants: list[Pattern],
        modifiers: RuleModifiers | None = None
    ) -> Generator[list[Pattern], Any, None]:
        pass


class FormBondTransform(ReactionTransformation):
    transformation = TransformationType.FORM_BOND
    bond: str
    patterns: list[Pattern]

    @classmethod
    def assert_correct_idxs(cls, center_in: list[tuple[int, int, int]],
                            center_out: list[tuple[int, int, int]]):
        # vreaking involves 2 different components
        assert len(center_in) == 2 and len(center_out) == 2

        center1_in, center2_in = center_in[0], center_in[1]
        idx1_in = center1_in[0]
        idx2_in = center2_in[0]

        center1_out, center2_out = center_out[0], center_out[1]
        idx1_out, mol_idx1_out, _ = center1_out
        idx2_out, mol_idx2_out, _ = center2_out

        # reaction centers in the input belong to 1 reagent
        # and in the output belong to 2 different reagents
        # and affected components in each center must be in different molecules
        assert idx1_in != idx2_in and idx1_out == idx2_out
        assert mol_idx1_out != mol_idx2_out

    def __init__(self,
                 reactants: list[Pattern],
                 graph_r: nx.Graph,
                 graph_p: nx.Graph,
                 center_in: list[tuple[int, int, int]],
                 center_out: list[tuple[int, int, int]]
                 ) -> None:
        super().__init__(graph_r)
        self.patterns = reactants

        FormBondTransform.assert_correct_idxs(center_in, center_out)

        center1_in, center2_in = center_in[0], center_in[1]
        center1_out, center2_out = center_out[0], center_out[1]

        idx1_in = center1_in[0]
        self.reactant_idx = idx1_in

        self.reaction_center_in = center_in
        self.reaction_center_out = center_out

        n1_in = graph_r.nodes[center1_in]
        n2_in = graph_r.nodes[center2_in]

        n1_out = graph_p.nodes[center1_out]
        n2_out = graph_p.nodes[center2_out]

        mol1_in_name = n1_in.get('molecule_name')
        mol1_out_name = n1_out.get('molecule_name')
        mol2_in_name = n2_in.get('molecule_name')
        mol2_out_name = n2_out.get('molecule_name')

        comp1_in_name = n1_in.get('comp_name', None)
        comp2_in_name = n2_in.get('comp_name', None)
        comp1_out_name = n1_out.get('comp_name', None)
        comp2_out_name = n2_out.get('comp_name', None)

        comp1_out_bond: str = n1_out.get('bond', '')
        comp2_out_bond: str = n2_out.get('bond', '')

        # assert molecule and component match between reaction centers in
        # input and output
        assert mol1_in_name == mol1_out_name
        assert mol2_in_name == mol2_out_name
        assert comp1_in_name == comp1_out_name
        assert comp2_in_name == comp2_out_name

        assert comp1_out_bond == comp2_out_bond and len(
            comp1_out_bond) and len(comp2_out_bond)

        self.bond = comp1_out_bond

    def apply(
        self,
        reactants: list[Pattern],
        modifiers: RuleModifiers | None = None
    ) -> Generator[list[Pattern], Any, None]:
        reactants = [reactant.copy() for reactant in reactants]

        center1 = self.reaction_center_in[0]
        center2 = self.reaction_center_in[1]
        i1 = center1[0]
        i2 = center2[0]

        specie1 = reactants[i1]
        specie2 = reactants[i2]

        match_func = sp_pattern_node_matching
        matcher1 = nx.isomorphism.GraphMatcher(
            specie1.graph, self.patterns[i1].graph, match_func)
        matcher2 = nx.isomorphism.GraphMatcher(
            specie2.graph, self.patterns[i2].graph, match_func)

        for map1 in matcher1.subgraph_isomorphisms_iter():
            map1 = {n2: n1 for n1, n2 in map1.items()}
            for map2 in matcher2.subgraph_isomorphisms_iter():
                map2 = {n2: n1 for n1, n2 in map2.items()}
                n1 = map1[center1[1:]]
                n2 = map2[center2[1:]]

                new_species = form_bond(specie1, specie2, n1, n2)

                products = reactants.copy()
                products[i1] = new_species
                products.pop(i2)

                yield products

    def __repr__(self) -> str:
        out = f"{self.transformation.value}: "
        return out


class BreakBondTransform(ReactionTransformation):
    transformation = TransformationType.BREAK_BOND
    patterns: list[Pattern]
    bond: str

    @classmethod
    def assert_correct_idxs(cls, center_in: list[tuple[int, int, int]],
                            center_out: list[tuple[int, int, int]]):
        # vreaking involves 2 different components
        assert len(center_in) == 2 and len(center_out) == 2

        center1_in, center2_in = center_in[0], center_in[1]
        idx1_in, mol_idx1_in, _ = center1_in
        idx2_in, mol_idx2_in, _ = center2_in

        center1_out, center2_out = center_out[0], center_out[1]
        idx1_out = center1_out[0]
        idx2_out = center2_out[0]

        # reaction centers in the input belong to 1 reagent
        # and in the output belong to 2 different reagents
        # and affected components in each center must be in different molecules
        assert idx1_in == idx2_in and idx1_out != idx2_out
        assert mol_idx1_in != mol_idx2_in

    def __init__(self,
                 reactants: list[Pattern],
                 graph_r: nx.Graph,
                 graph_p: nx.Graph,
                 center_in: list[tuple[int, int, int]],
                 center_out: list[tuple[int, int, int]]) -> None:
        super().__init__(graph_r)
        self.patterns = reactants

        BreakBondTransform.assert_correct_idxs(center_in, center_out)

        center1_in, center2_in = center_in[0], center_in[1]
        center1_out, center2_out = center_out[0], center_out[1]

        idx1_in = center1_in[0]
        self.reactant_idx = idx1_in

        self.reaction_center_in = center_in
        self.reaction_center_out = center_out

        n1_in = graph_r.nodes[center1_in]
        n2_in = graph_r.nodes[center2_in]

        n1_out = graph_p.nodes[center1_out]
        n2_out = graph_p.nodes[center2_out]

        mol1_in_name = n1_in.get('molecule_name')
        mol1_out_name = n1_out.get('molecule_name')
        mol2_in_name = n2_in.get('molecule_name')
        mol2_out_name = n2_out.get('molecule_name')

        comp1_in_name = n1_in.get('comp_name', None)
        comp2_in_name = n2_in.get('comp_name', None)
        comp1_out_name = n1_out.get('comp_name', None)
        comp2_out_name = n2_out.get('comp_name', None)

        comp1_in_bond: str = n1_in.get('bond', '')
        comp2_in_bond: str = n2_in.get('bond', '')

        # assert molecule and component match between reaction centers in
        # input and output
        assert mol1_in_name == mol1_out_name
        assert mol2_in_name == mol2_out_name
        assert comp1_in_name == comp1_out_name
        assert comp2_in_name == comp2_out_name

        assert comp1_in_bond == comp2_in_bond and len(
            comp1_in_bond) and len(comp2_in_bond)

        self.bond = comp1_in_bond

    def apply(
        self,
        reactants: list[Pattern],
        modifiers: RuleModifiers | None = None
    ) -> Generator[list[Pattern], Any, None]:
        reactants = [reactant.copy() for reactant in reactants]

        center1 = self.reaction_center_in[0]
        i = center1[0]
        specie = reactants[i]
        match_func = sp_pattern_node_matching
        matcher = nx.isomorphism.GraphMatcher(
            specie.graph, self.patterns_graph, match_func)
        
        for mapping in matcher.subgraph_isomorphisms_iter():
            mapping = {n2: n1 for n1, n2 in mapping.items()}
            n1 = mapping[center1]
            bond_id: str = specie.graph.nodes[n1]['bond']

            new_species = specie.break_bond(bond_id)
            products = reactants[:i] + list(new_species) + reactants[i+1:]

            yield products


    def __repr__(self) -> str:
        out = f"{self.transformation.value}:"
        return out


class ChangeStateAction(ReactionTransformation):
    transformation = TransformationType.CHANGE_COMPONENT_STATE
    # reaction_center_in: list[tuple[int, int, int]]
    # reaction_center_out: list[tuple[int, int, int]]
    patterns: list[Pattern]
    old_state: str
    new_state: str

    def __init__(self,
                 reactants: list[Pattern],
                 graph_r: nx.Graph,
                 graph_p: nx.Graph,
                 center_in: list[tuple[int, int, int]],
                 center_out: list[tuple[int, int, int]]) -> None:
        super().__init__(graph_r)
        self.patterns = reactants

        # only 1 component in the reaction center
        assert len(center_in) == 1 and len(center_out) == 1

        center1_in = center_in[0]
        center1_out = center_out[0]

        n1 = graph_r.nodes[center1_in]
        n2 = graph_p.nodes[center1_out]

        mol1_name = n1.get('molecule_name')
        comp1_name = n1.get('comp_name', None)
        comp1_state: str = n1.get('state', '')

        mol2_name = n2.get('molecule_name')
        comp2_name = n2.get('comp_name', None)
        comp2_state: str = n2.get('state', '')

        # assert valid state change
        assert mol1_name == mol2_name
        assert comp1_name == comp2_name
        assert (comp1_state != comp2_state
                and len(comp1_state) > 0
                and len(comp2_state) > 0)

        # reaction center
        self.reaction_center_in = center_in
        self.reaction_center_out = center_out

        # affected components
        self.old_state = comp1_state
        self.new_state = comp2_state

    def __repr__(self) -> str:
        comp1_str = f"{self.old_state}"
        comp2_str = f"{self.new_state}"
        out = (f"{self.transformation.value}: " +
               f"{self.reaction_center_in[0]} ({comp1_str})" +
               f" -> {self.reaction_center_out[0]} ({comp2_str})")
        return out

    def apply(
        self,
        reactants: list[Pattern],
        modifiers: RuleModifiers | None = None
    ) -> Generator[list[Pattern], Any, None]:
        # find affected reactant
        reactants = [reactant.copy() for reactant in reactants]
        # reactants_g = build_chemical_array_graph(products)
        # find map from reactant rules patterns graph to species reactants graph
        # note the reactant rules patterns graph is a isomorphic to a subgraph of
        # the species reactants graph
        center = self.reaction_center_in[0]
        i,j,k = center
        specie = reactants[i]

        match_func = sp_pattern_node_matching
        matcher = nx.isomorphism.GraphMatcher(
            specie.graph, self.patterns[i].graph, match_func)
        for mapping in matcher.subgraph_isomorphisms_iter():
            # reverse the map, we want the pattern to species map
            mapping = {n2: n1 for n1, n2 in mapping.items()}

            # affected node
            node: tuple[int, int] = mapping[self.reaction_center_in[0][1:]]
            j,k = node

            # change state and create new species
            molecules = [mol.copy() for mol in specie.molecules]
            molecules[j].change_state(k, self.new_state)
            new_specie = Pattern(molecules)

            products = reactants.copy()
            products[i] = new_specie
            yield products


def remove_dangling_bonds(mols: list[Molecule]):
    bonds: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)

    for i, mol in enumerate(mols):
        for j, comp in enumerate(mol.components):
            if comp.is_bonded() and not comp.is_bond_wildcard():
                bonds[comp.bond].append((i,j))

    for _, bonds_l in bonds.items():
        if len(bonds_l) == 2:
            continue

        # remove the dangling bonds
        for i, j in bonds_l:
            mols[i].components[j].bond = ''

    return mols


class DestroyMoleculeAction(ReactionTransformation):
    transformation = TransformationType.DESTROY_MOLECULE
    patterns: list[Pattern]
    react_idx: int

    def __init__(self,
                 reactants: list[Pattern],
                 graph_r: nx.Graph,
                 react_idx: int) -> None:
        super().__init__(graph_r)
        self.patterns = reactants
        self.react_idx = react_idx

    def __repr__(self) -> str:
        out = f"{self.transformation.value}: "
        return out

    def apply(
        self,
        reactants: list[Pattern],
        modifiers: RuleModifiers | None = None
    ) -> Generator[list[Pattern], Any, None]:
        # find affected reactant
        reactants = [reactant.copy() for reactant in reactants]

        idx = self.react_idx
        pattern = self.patterns[idx]
        if not modifiers or not modifiers.delete_molecules:
            products = reactants.copy()
            products.pop(idx)
            yield products
        else:
            affected_specie = reactants[idx].copy()

            # we need to find the isomorphism between the affected reactant and the pattern
            # and then delete the correct molecules and any dangling bonds
            match_func = sp_pattern_node_matching
            matcher = nx.isomorphism.GraphMatcher(
                affected_specie.graph, pattern.graph, match_func)
            for mapping in matcher.subgraph_isomorphisms_iter():
                products = reactants.copy()

                # remove affected molecules
                n_mols = affected_specie.number_of_molecules()
                idxs_to_remove = set(n1[0] for n1 in mapping.keys())
                mols = [affected_specie.molecules[i].copy() for i in range(n_mols)
                        if i not in idxs_to_remove]

                # remove dangling bonds
                mols = remove_dangling_bonds(mols)

                if len(mols) > 0:
                    new_pattern = Pattern(mols)
                    # if deleting the molecules results into multiple disconnected molecules
                    if not new_pattern.is_connected():
                        new_patterns = split_pattern_into_connected_components(new_pattern)
                        products = products[:idx] + new_patterns + products[idx+1:]
                    else:
                        products[idx] = new_pattern
                else:
                    products.pop(idx)
                yield products
                return


class CreateMoleculeAction(ReactionTransformation):
    transformation = TransformationType.CREATE_MOLECULE
    patterns: list[Pattern]
    created_pattern: Pattern

    def __init__(self,
                 reactants: list[Pattern],
                 graph_r: nx.Graph,
                 created_pattern: Pattern) -> None:
        super().__init__(graph_r)
        self.patterns = reactants
        self.created_pattern = created_pattern

    def __repr__(self) -> str:
        out = f"{self.transformation.value}: "
        return out

    def apply(
        self,
        reactants: list[Pattern],
        modifiers: RuleModifiers | None = None
    ) -> Generator[list[Pattern], Any, None]:
        # find affected reactant
        reactants = [reactant.copy() for reactant in reactants]

        products = reactants.copy()
        products.append(self.created_pattern.copy())
        yield products


def apply_transforms(
    reactants: list[Pattern],
    actions: list[ReactionTransformation],
    modifiers: RuleModifiers | None = None
) -> Generator[list[Pattern], Any, None]:
    # Start with the initial list of reactants as the first generator
    def generator_chain(reactants: list[Pattern], actions: list[ReactionTransformation]) -> Generator[list[Pattern], Any, None]:
        # Initial generator yields the initial reactants
        res = [reactants]

        # Apply each action to the results of the previous action
        for action in actions:
            # Apply the current action to each result from the previous step
            res = [p for r in res for p in action.apply(r, modifiers)]

        # Yield the final results after all actions have been applied
        yield from res

    return generator_chain(reactants, actions)
