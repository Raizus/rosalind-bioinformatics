import pytest

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Observable import Observable
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern


class TestObservable():
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

    @pytest.mark.parametrize("declaration, expected", [
        ("Molecules A A(b)", ('Molecules', 'A')),
        ("Species B B(a)", ('Species', 'B')),
        ("Molecules C C(a!?)", ('Molecules', 'C')),
        ("Species test_2 A(b!+)", ('Species', 'test_2')),
    ])
    def test_valid(self, declaration: str, expected: tuple[str, str]):
        observable = Observable.from_declaration(declaration, self.molecules)
        assert observable.type == expected[0]
        assert observable.label == expected[1]

    @pytest.mark.parametrize("declaration, expected_num_matches, total_counts", [
        ("Molecules A A(b)", 4, 6),
        ("Molecules A A(b,b)", 2, 2),
        ("Species A A(b)", 4, 4),
        ("Molecules A_bonded A(b!+)", 4, 6),
        ("Species A_bonded A(b!+)", 4, 4),
        ("Molecules C_bonded C(a!+)", 3, 3),
        ("Molecules C C(a)", 1, 1),
    ])
    def test_pattern_matching(self, declaration: str, expected_num_matches: int, total_counts: int):
        observable = Observable.from_declaration(declaration, self.molecules)
        species = [Pattern.from_declaration(
            species_str, self.molecules) for species_str in self.species_str_list]
        matches = {pattern2: observable.match_species(
            pattern2) for pattern2 in species}
        filtered_matches = {patt: count for patt, count in matches.items() if count > 0}
        assert len(filtered_matches) == expected_num_matches
        assert sum(filtered_matches.values()) == total_counts
