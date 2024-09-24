import abc
from enum import Enum
from functools import reduce

from BioInfoToolkit.RuleBasedModel.model.Component import Component
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern


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


class ReactionTransformation(abc.ABC):
    transformation: TransformationType
    # tuple with (reactant_idx, molecule_idx, component_idx)
    # list of components that participate in a tranformation
    #   - Form bond: involves 2 components in two different molecules and different reactants
    #   - Break bond: involves 2 components in two different molecules in the same reactants
    #   - Change Component State: Involves a single component

    reaction_center_in: list[tuple[int, int, int]]
    reaction_center_out: list[tuple[int, int, int]]

    # list of components that do not participate in the transformation,
    # but are nonetheless required and influence the reaction rate
    reaction_context_in: list[tuple[int, int, int]]
    reaction_context_out: list[tuple[int, int, int]]

    def __init__(self, reactants: list[Pattern], products: list[Pattern]) -> None:
        self.reaction_center_in = []
        self.reaction_center_out = []
        self.reaction_context_in = idxs_from_reactants(reactants)
        self.reaction_context_out = idxs_from_reactants(products)

    @abc.abstractmethod
    def __call__(self, reactants: list[Pattern]) -> list[Pattern]:
        pass


class FormBondTransform(ReactionTransformation):
    transformation = TransformationType.FORM_BOND
    bond: str

    def __init__(self,
                 reactants: list[Pattern],
                 products: list[Pattern],
                 center_in: list[tuple[int, int, int]],
                 center_out: list[tuple[int, int, int]]
                 ) -> None:

        assert len(center_in) == 2 and len(center_out) == 2

        reactants_idxs = idxs_from_reactants(reactants)
        product_idxs = idxs_from_reactants(products)

        center1_in, center2_in = center_in[0], center_in[1]
        idx1_in, mol_idx1_in, comp_idx1_in = center1_in
        idx2_in, mol_idx2_in, comp_idx2_in = center2_in

        center1_out, center2_out = center_out[0], center_out[1]
        idx1_out, mol_idx1_out, comp_idx1_out = center1_out
        idx2_out, mol_idx2_out, comp_idx2_out = center2_out

        assert idx1_in != idx2_in and idx1_out == idx2_out
        assert mol_idx1_in != mol_idx2_in and mol_idx1_out != mol_idx2_out

        self.reaction_center_in = center_in
        self.reaction_center_out = center_out

        for idx in center_in:
            reactants_idxs.remove(idx)
        for idx in center_out:
            product_idxs.remove(idx)

        self.reaction_context_in = reactants_idxs
        self.reaction_context_out = product_idxs

        # find bond label
        reactant1_in, molecule1_in, comp1_in = objs_from_idx(
            reactants, center1_in)
        reactant2_in, molecule2_in, comp2_in = objs_from_idx(
            reactants, center2_in)
        reactant1_out, molecule1_out, comp1_out = objs_from_idx(
            products, center1_out)
        reactant2_out, molecule2_out, comp2_out = objs_from_idx(
            products, center2_out)

        assert molecule1_in._name == molecule1_out._name
        assert molecule2_in._name == molecule2_out._name
        assert comp1_in.name == comp1_out.name
        assert comp2_in.name == comp2_out.name

        bond1 = comp1_out.bond
        bond2 = comp2_out.bond

        assert bond1 == bond2 and len(bond1) and len(bond2)
        self.bond = bond1

    # def apply(self, reactants: list[Pattern]) -> list[Pattern]:
    #     reaction_center_in = self.reaction_center_in
    #     reaction_center_out = self.reaction_center_out
    
    #     center1_in, center2_in = reaction_center_in[0], reaction_center_in[1]
    #     reactant1_in, molecule1_in, comp1_in = objs_from_idx(
    #         reactants, center1_in)
    #     reactant2_in, molecule2_in, comp2_in = objs_from_idx(
    #         reactants, center2_in)

    #     # need to bond reactant1 to reactant2 and set the bond in comp1 and comp2


class ChangeStateAction(ReactionTransformation):
    transformation = TransformationType.CHANGE_COMPONENT_STATE
    # reaction_center_in: list[tuple[int, int, int]]
    # reaction_center_out: list[tuple[int, int, int]]
    component_in: Component
    component_out: Component

    def __init__(self,
                 reactants: list[Pattern],
                 products: list[Pattern],
                 center_in: list[tuple[int, int, int]],
                 center_out: list[tuple[int, int, int]]) -> None:
        super().__init__(reactants, products)

        # only 1 component in the reaction center
        assert len(center_in) == 1 and len(center_out) == 1

        center1_in = center_in[0]
        center1_out = center_out[0]

        _, molecule1_in, comp1_in = objs_from_idx(
            reactants, center1_in)
        _, molecule1_out, comp1_out = objs_from_idx(
            products, center1_out)

        # assert valid state change
        assert molecule1_in._name == molecule1_out._name
        assert comp1_in.name == comp1_out.name
        assert (comp1_in.state != comp1_out.state
                and len(comp1_in.state) > 0 
                and len(comp1_out.state) > 0)

        # center and context
        self.reaction_center_in = center_in
        self.reaction_center_out = center_out
        self.reaction_context_in = [idx for idx in self.reaction_context_in if idx not in center_in]
        self.reaction_context_out = [
            idx for idx in self.reaction_context_out if idx not in center_out]

        # affected components
        self.component_in = comp1_in
        self.component_out = comp1_out

    def __repr__(self) -> str:
        out = (f"{self.transformation.value}: {self.reaction_center_in[0]} ({self.component_in})" +
            f" -> {self.reaction_center_out[0]} ({self.component_out})")
        return out

    def __call__(self, reactants: list[Pattern]) -> list[Pattern]:
        # find affected reactant
        products = [reactant.copy() for reactant in reactants]
        i,j,k = self.reaction_center_in[0]
        products[i].update_component_state((j, k), self.component_out.state)

        return products


class CreateMoleculeAction(ReactionTransformation):
    transformation = TransformationType.CREATE_MOLECULE
    molecule_pattern: Pattern
    product_idx: int

    def __init__(self, reactants: list[Pattern],
                 products: list[Pattern],
                 center_out: list[tuple[int, int, int]]) -> None:
        super().__init__(reactants, products)

        assert len(center_out) > 0
        assert all(center_out[0][0] == idx[0] for idx in center_out)
        assert all(center_out[0][1] == idx[1] for idx in center_out)
        i, j = center_out[0][0], center_out[0][1]
        self.reaction_center_out = center_out

        # verify center (must be all components belonging to the same molecule and the same pattern)
        extracted_context = [idx for idx in self.reaction_context_out if idx[0] == i]
        assert set(center_out) == set(extracted_context)

        # remove center from context
        self.reaction_context_out = [idx for idx in self.reaction_context_out
                                     if idx not in extracted_context]

        # extract molecule from center
        self.molecule_pattern = products[i].copy()
        self.product_idx = i

    def __call__(self, reactants: list[Pattern]) -> list[Pattern]:
        products = [reactant.copy() for reactant in reactants]
        i = self.product_idx

        # make sure we insert in the correct index
        products = products[:i] + [self.molecule_pattern] + products[i+1:]
        return products


def apply_transforms(reactants: list[Pattern], actions: list[ReactionTransformation]) -> list[Pattern]:
    products = reduce(lambda x, y: y(x), actions, reactants)
    return products
