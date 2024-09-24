import pytest

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Pattern import generate_species


def test_valid_1():
    declaration = "A(x~a~b,x~a~b,y)"
    molecule = MoleculeType(declaration)

    assert isinstance(molecule, MoleculeType)


def test_invalid_1():
    declaration = "A(x~a~b,x~b~c,y)"
    with pytest.raises(ValueError):
        MoleculeType(declaration)


def test_generate_species_1():
    declaration = "A(x~a~b,x~a~b,y)"
    molecule = MoleculeType(declaration)

    species = [x for x in generate_species(molecule)]

    assert len(species) == 4