import pytest

from BioInfoToolkit.RuleBasedModel.model.Component import Component
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Pattern import Molecule, Pattern, match_patterns


def test_valid_1():
    components = [Component('x', {'a', 'b'}),
                  Component('x', {'a', 'b'}),
                  Component('y', set())]
    reactant = Molecule('A', components)

    assert isinstance(reactant, Molecule)


def test_from_declaration():
    molecules: dict[str, MoleculeType] = {
        "A": MoleculeType("A(x~a~b,x~a~b,y)")
    }

    declaration = "A(x~a,x~a,y)"
    reactant = Molecule.from_declaration(declaration, molecules)
    assert isinstance(reactant, Molecule)


class TestMolecule():
    molecules: dict[str, MoleculeType] = {
        "A": MoleculeType("A(x~a~b,x~a~b,y)")
    }

    @pytest.mark.parametrize("declaration", [
        "A()", 
        "A(x~a,x~a,y)", "A(x~a,x~b,y)", "A(x~b,x~a,y)", "A(x~b,x~b,y)",
        "A(y,x~a,x~a)", "A(x~a,y,x~b)", "A(x~b,y,x~a)", "A(y,x~b,x~b)",
        "A(x~a,x~a)", "A(x~a,y)", "A(x~b,y)", "A(y)"
    ])
    def test_valid(self, declaration: str):
        reactant = Molecule.from_declaration(declaration, self.molecules)
        assert isinstance(reactant, Molecule)

    @pytest.mark.parametrize("declaration", [
        "A(z)",
        "A(x~c)", "A(y~a)", 
        "A(x~a,x~a,x~b)", "A(x~a, y, y)",
        "B()"
    ])
    def test_invalid(self, declaration: str):
        with pytest.raises(ValueError):
            Molecule.from_declaration(declaration, self.molecules)

    @pytest.mark.parametrize("declaration", [
        "A(x~a,x~a,y)", "A(x~a,x~b,y)", "A(x~b,x~a,y)", "A(x~b,x~b,y)",
        "A(y,x~a,x~a)", "A(x~a,y,x~b)", "A(x~b,y,x~a)", "A(y,x~b,x~b)",
        "A(x~a,x~a)"
    ])
    def test_is_species(self, declaration: str):
        reactant = Molecule.from_declaration(declaration, self.molecules)
        assert reactant.is_specie() is True

    @pytest.mark.parametrize("declaration", [
        "A()",
        "A(x~a,y)", "A(x~b,y)", "A(y)"
    ])
    def test_is_pattern(self, declaration: str):
        molecule = Molecule.from_declaration(declaration, self.molecules)
        assert molecule.is_pattern() is True

    @pytest.mark.parametrize("declaration, expected", [
        ("A()", 4),
        ("A(y)", 4),
        ("A(x~a,y)", 2), 
        ("A(x~b,y)", 2), 
        ("A(x~b)", 2),
        ("A(x~a,x~b,y)", 1),
        ("A(x~a,x~a,y)", 1),
        ("A(x~a,x~a)", 1),
    ])
    def test_generate_species(self, declaration: str, expected: int):
        reactant = Molecule.from_declaration(declaration, self.molecules)
        species = [specie for specie in reactant.generate_species()]
        assert len(species) == expected


class TestPattern():
    molecules: dict[str, MoleculeType] = {
        "L": MoleculeType("L(rec)"),
        "R": MoleculeType("R(lig,ch~closed~open)")
    }

    molecules2: dict[str, MoleculeType] = {
        "A": MoleculeType("A(x,y)"),
        "B": MoleculeType("B(p,q~a~b)"),
        "C": MoleculeType("C(r~d~e,s)"),
        "D": MoleculeType("D(a,a)")
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
        complex = Pattern.from_declaration(declaration, self.molecules2)
        assert isinstance(complex, Pattern)

    # @pytest.mark.parametrize("declaration", [
    #     "L(rec!0).R(lig!0,ch~open)", "L(rec!0).R(lig!0,ch~closed)"
    # ])
    # def test_is_species(self, declaration: str):
    #     complex = ComplexReactant.from_declaration(declaration, self.molecules)
    #     assert complex.is_specie() == True

    # @pytest.mark.parametrize("declaration", [
    #     "L(rec!0).R(lig!0)", "R(lig!0).L(rec!0)", 
    #     "L(rec!1).R(lig!1)", "R(lig!1).L(rec!1)"
    # ])
    # def test_is_pattern(self, declaration: str):
    #     complex = ComplexReactant.from_declaration(declaration, self.molecules)
    #     assert complex.is_pattern() == True

    # @pytest.mark.parametrize("declaration, expected", [
    #     ("L(rec!0).R(lig!0)", 2),
    #     ("R(lig!1).L(rec!1)", 2),
    #     ("L(rec!0).R(lig!0,ch~open)", 1),
    #     ("L(rec!0).R(lig!0,ch~closed)", 1),
    # ])
    # def test_generate_species(self, declaration: str, expected: int):
    #     complex = ComplexReactant.from_declaration(declaration, self.molecules)
    #     species = [specie for specie in complex.generate_species()]
    #     assert len(species) == expected

    # @pytest.mark.parametrize("declaration, expected", [
    #     ("D(a!0).A(x!0)", 2),
    # ])
    # def test_generate_species_2(self, declaration: str, expected: int):
    #     complex = ComplexReactant.from_declaration(declaration, self.molecules2)
    #     species = [specie for specie in complex.generate_species()]
    #     assert len(species) == expected

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

        match_counts = match_patterns(pattern1, pattern2)
        match = match_counts > 0
        assert match == expected
