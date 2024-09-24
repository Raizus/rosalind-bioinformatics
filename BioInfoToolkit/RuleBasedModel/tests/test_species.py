import pytest

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern, match_patterns
from BioInfoToolkit.RuleBasedModel.model.Species import Species, generate_species_from_molecule_type


class TestSpecies():
    molecules: dict[str, MoleculeType] = {
        "A": MoleculeType("A(b,b,c)"),
        "B": MoleculeType("B(a)"),
        "C": MoleculeType("C(a)"),
    }

    species_str_list = [
        "A(b,b,c)", "B(a)", "C(a)",
        # "A(b,c,b)",
        "A(b!0,b,c).B(a!0)",
        "A(b!0,b!1,c).B(a!0).B(a!1)",
        "A(b,b,c!2).C(a!2)",
        "A(b!0,b,c!2).B(a!0).C(a!2)",
        "A(b!0,b!1,c!2).B(a!0).B(a!1).C(a!2)"
    ]

    @pytest.mark.parametrize("declaration", species_str_list)
    def test_valid(self, declaration: str):
        specie = Pattern.from_declaration(declaration, self.molecules)
        assert specie is not None

    @pytest.mark.parametrize("declaration, count", [
        ('A(b)', 4),
        ('A(b!+)', 4),
        ('A(b!?)', 6),
        ('A(b).C()', 2)
    ])
    def test_matching(self, declaration: str, count: int):
        target_pattern = Pattern.from_declaration(declaration, self.molecules)
        patterns = [Pattern.from_declaration(decl, self.molecules) for decl in self.species_str_list]
        matches = [pattern2 for pattern2 in patterns if match_patterns(target_pattern, pattern2)]
        l = len(matches)
        assert l == count


# class TestSpecies2():
#     molecules = {
#         "L": MoleculeType("L(t)"),
#         "T": MoleculeType("T(l,r,Phos~U~P,Meth~A~B~C)"),
#         "CheY": MoleculeType("CheY(Phos~U~P)"),
#         "CheZ": MoleculeType("CheZ()"),
#         "CheB": MoleculeType("CheB(Phos~U~P)"),
#         "CheR": MoleculeType("CheR(t)")
#     }

#     species_str_list = [
#         "L(t)",
#         "T(Phos~U, l)",
#         "T(Phos~P, l)",
#         "CheY(Phos~U)",
#         "CheY(Phos~P)",
#         "CheZ()",
#         "L(t!1).T(Phos~U, l!1)",
#         "L(t!1).T(Phos~P, l!1)",
#     ]

#     @pytest.mark.parametrize("declaration", species_str_list)
#     def test_valid(self, declaration: str):
#         specie = Pattern.from_declaration(declaration, self.molecules)
#         assert specie is not None

#     @pytest.mark.parametrize("declaration, count", [
#         ('A(b)', 4),
#         ('A(b!+)', 4),
#         ('A(b!?)', 6),
#         ('A(b).C()', 2)
#     ])
#     def test_matching(self, declaration: str, count: int):
#         target_pattern = Pattern.from_declaration(declaration, self.molecules)
#         patterns = [Pattern.from_declaration(
#             decl, self.molecules) for decl in self.species_str_list]
#         matches = [pattern2 for pattern2 in patterns if match_patterns(
#             target_pattern, pattern2)]
#         l = len(matches)
#         assert l == count



class TestSpeciesGeneration():
    molecules: dict[str, MoleculeType] = {
        "L": MoleculeType("L(t)"),
        "T": MoleculeType("T(l,r,Phos~U~P,Meth~A~B~C)"),
        "CheY": MoleculeType("CheY(Phos~U~P)"),
        "CheZ": MoleculeType("CheZ()"),
        "CheB": MoleculeType("CheB(Phos~U~P)"),
        "CheR": MoleculeType("CheR(t)"),
    }

    def test_species_generation(self):
        species = [specie 
                   for molecule_type in self.molecules.values()
                   for specie in generate_species_from_molecule_type(molecule_type)
                   ]
        assert len(species) == 13
