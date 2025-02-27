"""
    Tests MoleculeType declarations
"""
import pytest
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType


class TestMoleculeType:
    """
        Tests for valid and invalid MoleculeType declarations
    """
    @pytest.mark.parametrize("declaration", [
        "A()",
        "A( )",
        "A(y)",
        "A(y, y, y)",
        "A(y,y,y)",
        "A(x~a~b,x~a~b,y)",
        "Abc()",
    ])
    def test_valid(self, declaration: str):
        """Tests valid MoleculeType declarations

        Args:
            declaration (str): _description_
        """
        molecule = MoleculeType.from_declaration(declaration)
        assert isinstance(molecule, MoleculeType)

    @pytest.mark.parametrize("declaration", [
        "A(x~a~b,x~b~c,y)",
    ])
    def test_invalid(self, declaration: str):
        """Tests for invalid MoleculeType declarations

        Args:
            declaration (str): _description_
        """
        declaration = "A(x~a~b,x~b~c,y)"
        with pytest.raises(ValueError):
            MoleculeType.from_declaration(declaration)
