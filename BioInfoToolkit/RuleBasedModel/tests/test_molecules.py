import pytest

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType


class TestMoleculeType:

    @pytest.mark.parametrize("declaration", [
        "A(x~a~b,x~a~b,y)",
        "A(y)",
        "Abc()",
    ])
    def test_valid(self, declaration: str):
        molecule = MoleculeType.from_declaration(declaration)
        assert isinstance(molecule, MoleculeType)

    @pytest.mark.parametrize("declaration", [
        "A(x~a~b,x~b~c,y)",
    ])
    def test_invalid(self, declaration: str):
        declaration = "A(x~a~b,x~b~c,y)"
        with pytest.raises(ValueError):
            MoleculeType.from_declaration(declaration)
