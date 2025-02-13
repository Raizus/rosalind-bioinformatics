import pytest

from BioInfoToolkit.RuleBasedModel.model.Component import Component
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Pattern import Molecule, Pattern, match_pattern_specie


def test_valid_molecule():
    components = [Component('x', {'a', 'b'}),
                  Component('x', {'a', 'b'}),
                  Component('y', set())]
    reactant = Molecule('A', components)
    assert isinstance(reactant, Molecule)


def test_from_declaration():
    molecules: dict[str, MoleculeType] = {
        "A": MoleculeType.from_declaration("A(x~a~b,x~a~b,y)")
    }

    declaration = "A(x~a,x~a,y)"
    reactant = Molecule.from_declaration(declaration, molecules)
    assert isinstance(reactant, Molecule)


class TestMolecule():
    molecule_types: dict[str, MoleculeType] = {
        "A": MoleculeType.from_declaration("A(x~a~b,x~a~b,y)")
    }

    @pytest.mark.parametrize("declaration", [
        "A()", 
        "A(x~a,x~a,y)", "A(x~a,x~b,y)", "A(x~b,x~a,y)", "A(x~b,x~b,y)",
        "A(y,x~a,x~a)", "A(x~a,y,x~b)", "A(x~b,y,x~a)", "A(y,x~b,x~b)",
        "A(x~a,x~a)", "A(x~a,y)", "A(x~b,y)", "A(y)"
    ])
    def test_valid(self, declaration: str):
        reactant = Molecule.from_declaration(declaration, self.molecule_types)
        assert isinstance(reactant, Molecule)

    @pytest.mark.parametrize("declaration", [
        "A(z)",
        "A(x~c)", "A(y~a)", 
        "A(x~a,x~a,x~b)",
        "B()",
        # "A(x~a, y, y)",
    ])
    def test_invalid(self, declaration: str):
        with pytest.raises(ValueError):
            Molecule.from_declaration(declaration, self.molecule_types)


class TestPattern():
    molecules: dict[str, MoleculeType] = {
        "L": MoleculeType.from_declaration("L(rec)"),
        "R": MoleculeType.from_declaration("R(lig,ch~closed~open)")
    }

    molecules2: dict[str, MoleculeType] = {
        "A": MoleculeType.from_declaration("A(x,y)"),
        "B": MoleculeType.from_declaration("B(p,q~a~b)"),
        "C": MoleculeType.from_declaration("C(r~d~e,s)"),
        "D": MoleculeType.from_declaration("D(a,a)")
    }

    @pytest.mark.parametrize("declaration", [
        "L(rec!0).R(lig!0)", "R(lig!0).L(rec!0)", "L(rec!1).R(lig!1)", "R(lig!1).L(rec!1)",
        "L(rec!0).R(lig!0,ch~open)", "L(rec!0).R(lig!0,ch~closed)",
        "R(lig!0,ch~open).L(rec!0)", "R(lig!0,ch~closed).L(rec!0)",
    ])
    def test_valid_1(self, declaration: str):
        pattern = Pattern.from_declaration(declaration, self.molecules)
        assert isinstance(pattern, Pattern)

    @pytest.mark.parametrize("declaration", [
        "A(x,y!0).B(p!0,q~a!1).C(r~d!1,s!2).A(x!2,y)",
        "B(p!0,q~a!1).C(r~d!1,s!2).A(x!2,y).A(x,y!0)"
    ])
    def test_valid_2(self, declaration: str):
        pattern = Pattern.from_declaration(declaration, self.molecules2)
        assert isinstance(pattern, Pattern)

    @pytest.mark.parametrize("pattern1_str, pattern2_str, expected", [
        ("R(lig!0,ch~open)", "R(lig!0,ch~open)", True),     # exactly equal
        ("R(ch~open,lig!0)", "R(lig!0,ch~open)", True),     # permute components
        ("R(lig!0,ch~open)", "R(lig!?,ch~open)", True),     # bond wildmark '?'
        ("R(lig!?,ch~open)", "R(lig!+,ch~open)", True),     # bond wildmark '?' and '+'
        ("R(lig!0,ch~open)", "R(lig!+,ch~open)", True),     # bond wildmark '+'
        ("R(lig!0,ch~open)", "R(lig!0,ch~closed)", False),  # state change
        ("R(lig,ch)", "R(lig!0,ch)", False),                # with vs without bond

        # omited vs explicit component, component mismatch
        ("R(lig!0,ch)", "R(lig!0)", False),
        ("R(lig,ch~open)", "R(lig!0,ch~open)", False),
        ("R(lig,ch)", "R(lig)", False),

        ("L(rec!0).R(lig!0)", "L(rec!0).R(lig!0)", True),   # exact equality
        ("L(rec!0).R(lig!0)", "R(lig!0).L(rec!0)", True),   # permute molecules
        ("L(rec!0).R(lig!0)", "R(lig!2).L(rec!2)", True),   # check different bond labels

        # different bond labels but missing edge in the graph, so does not match)
        ("L(rec!?).R(lig!0)", "R(lig!0).L(rec!0)", False),
        ("L(rec!+).R(lig!0)", "R(lig!0).L(rec!0)", False),
    ])
    def test_equality(self, pattern1_str: str, pattern2_str: str, expected: bool):
        pattern1 = Pattern.from_declaration(pattern1_str, self.molecules)
        pattern2 = Pattern.from_declaration(pattern2_str, self.molecules)
        match = pattern1 == pattern2
        assert match == expected

    @pytest.mark.parametrize("pattern1_str, pattern2_str, expected", [
        ("R(lig!0).L(rec!0)", "R(lig!0,ch~open).L(rec!0)", True),
        ("R(lig!0,ch~open).L(rec!0)", "R(lig!0).L(rec!0)", False)
    ])
    def test_pattern_matching(self, pattern1_str: str, pattern2_str: str, expected: bool):
        pattern1 = Pattern.from_declaration(pattern1_str, self.molecules)
        pattern2 = Pattern.from_declaration(pattern2_str, self.molecules)

        match_counts = match_pattern_specie(pattern1, pattern2)
        match = match_counts > 0
        assert match == expected
