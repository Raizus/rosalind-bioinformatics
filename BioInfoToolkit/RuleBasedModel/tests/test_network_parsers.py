import pytest

from BioInfoToolkit.RuleBasedModel.network.parsers import parse_observables_group, parse_parameters, parse_reaction, parse_seed_species


class TestParseParameters():
    @pytest.mark.parametrize("declaration, p_id, label, expression", [
        ("1 NaV2                 6.02e8  # Constant", 1, 'NaV2', "6.02e8"),
        ("6 k_lr_bind            8.8e6/NaV2", 6, 'k_lr_bind', "8.8e6/NaV2"),
        ("6 Lig_tot  lig_conc*Na*Vec  # units: molecules",
         6, 'Lig_tot', 'lig_conc*Na*Vec')
    ])
    def test_valid(self, declaration: str, p_id: int, label: str, expression: str):
        parsed = parse_parameters(declaration)

        assert parsed['id'] == p_id
        assert parsed["name"] == label
        assert parsed['expression'] == expression


class TestParseObservablesGroup():
    @pytest.mark.parametrize("declaration, g_id, label, species", [
        ("1 LynFree              2",
         1, 'LynFree', [(2,1)]),
        ("2 RecMon               4,5,6,7,13,14,15,338",
         2, 'RecMon', [(4, 1), (5, 1), (6, 1), (7, 1), (13, 1), (14, 1), (15, 1), (338, 1)]),
        ("4 RecPbeta             18,20,21,2*23,42,44,2*45,2*46,47",
         4, 'RecPbeta',
         [(18, 1), (20, 1), (21, 1), (23, 2), (42, 1), (44, 1), (45, 2), (46, 2), (47, 1)]),
    ])
    def test_valid(self, declaration: str, g_id: int, label: str, species: list[tuple[int,int]]):
        parsed = parse_observables_group(declaration)

        assert parsed['id'] == g_id
        assert parsed["label"] == label
        assert parsed['species'] == species


class TestParseReactions():
    @pytest.mark.parametrize("declaration, r_id, reactants, products, rate", [
        ("1 1,4 5 2*kp1  # R1", 1, [1, 4], [5], "2*kp1"),
        ("2 2,4 6 kpL  # R3", 2, [2, 4], [6], "kpL"),
        ("3 1,6 7 2*kp1  # R1", 3, [1, 6], [7], "2*kp1"),
        ("4 5 1,4 km1  # _reverse_R1", 4, [5], [1, 4], "km1"),
        ("5 4,5 8 kp2  # R2", 5, [4, 5], [8], "kp2"),
    ])
    def test_valid(self, declaration: str, r_id: int, 
                   reactants: list[int], products: list[int],
                   rate: str):
        parsed = parse_reaction(declaration)

        assert parsed['id'] == r_id
        assert parsed["reactants"] == reactants
        assert parsed["products"] == products
        assert parsed["rate"] == rate


class TestSeedSpeciesParser():
    @pytest.mark.parametrize("declaration, sp_id, expression", [
        ("1 L(t) L0", 1, "L0"),
        ("2 T(Phos~U,l) _InitialConc1", 2, "_InitialConc1"),
        ("3 T(Phos~P,l) _InitialConc2", 3, "_InitialConc2"),
        ("4 CheY(Phos~U) _InitialConc3", 4, "_InitialConc3"),
        ("5 CheY(Phos~P) _InitialConc4", 5, "_InitialConc4"),
        ("6 CheZ() CheZ0", 6, "CheZ0"),
        ("7 L(t!1).T(Phos~U,l!1) 0", 7, "0"),
        ("8 L(t!1).T(Phos~P,l!1) 0", 8, "0"),
    ])
    def test_valid(self, declaration: str, sp_id: int, expression: str):
        parsed = parse_seed_species(declaration)

        assert parsed['id'] == sp_id
        assert parsed["expression"] == expression
