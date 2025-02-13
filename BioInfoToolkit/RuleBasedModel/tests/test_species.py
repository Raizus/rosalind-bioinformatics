import pytest

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern, match_pattern_specie


class TestSpecies():
    molecules: dict[str, MoleculeType] = {
        "A": MoleculeType.from_declaration("A(b,b,c)"),
        "B": MoleculeType.from_declaration("B(a)"),
        "C": MoleculeType.from_declaration("C(a)"),
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
        patterns = [Pattern.from_declaration(decl, self.molecules)
                    for decl in self.species_str_list]
        matches = [pattern2 for pattern2 in patterns
                   if match_pattern_specie(target_pattern, pattern2)]
        l = len(matches)
        assert l == count


class TestPatternMatching:
    molecules: dict[str, MoleculeType] = {
        "L": MoleculeType.from_declaration("L(l,l)"),
        "R": MoleculeType.from_declaration("R(r,r)"),
    }

    @pytest.mark.parametrize("pattern_str, species_str", [
        ('R(r!1).L(l,l!1)', 'L(l!1,l!2).L(l,l!3).R(r!1,r!3).R(r,r!2)'),
    ])
    def test_matching(self, pattern_str: str, specie_str: str):
        pattern = Pattern.from_declaration(pattern_str, self.molecules)
        specie = Pattern.from_declaration(specie_str, self.molecules)

        assert match_pattern_specie(pattern, specie) == 1
