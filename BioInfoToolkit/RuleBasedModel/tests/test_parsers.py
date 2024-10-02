import pytest

from BioInfoToolkit.RuleBasedModel.utils.model_parsers import parse_observable, \
    parse_pattern, parse_reactants_sum, parse_reaction_rule, parse_molecule, parse_molecule_type


class TestParseMoleculeType():
    @pytest.mark.parametrize("declaration, name, count", [
        ("A()", "A", 0),
        ("ABC()", "ABC", 0),
        ("A(x)", "A", 1),
        ("A(x,x,y)", "A", 3),
        ("A(x~a~b,x~a~b,y)", "A", 3),
        ("A(x~a~b~c,x~a~b~c,y)", "A", 3),
        (" A(x~a~b~c,x~a~b~c,y) ", "A", 3)
    ])
    def test_valid(self, declaration: str, name: str, count: int):
        parsed = parse_molecule_type(declaration)
        assert parsed['name'] == name
        assert len(parsed["components"]) == count

    def test_states(self):
        declaration = "A(x~a~b~c,x~a~b~c,y)"
        parsed = parse_molecule_type(declaration)
        assert parsed['name'] == "A"
        assert len(parsed["components"]) == 3

        assert parsed["components"][0]["name"] == 'x'
        assert parsed["components"][0]["states"] == set(['a', 'b', 'c'])
        assert parsed["components"][1]["name"] == 'x'
        assert parsed["components"][1]["states"] == set(['a', 'b', 'c'])
        assert parsed["components"][2]["name"] == 'y'
        assert parsed["components"][2]["states"] == set()

    @pytest.mark.parametrize("declaration", [
        "A( )", "ABC( )",
        "A(x~a~b, x~a~b,y)", "C( x)", "C2(y )"
    ])
    def test_invalid(self, declaration: str):
        with pytest.raises(ValueError):
            parse_molecule_type(declaration)


class TestParseMolecule():
    @pytest.mark.parametrize("declaration", [
        "A()", "ABC()", "A() ",
        "A(x~a,x~a,y)", "A(x1~a!0,x1~b,y!1)", "ABC(x1,x2~abc,x3!0,x4~a!1)",
    ])
    def test_valid(self, declaration: str):
        parsed = parse_molecule(declaration)
        assert parsed is not None

    def test_components(self):
        declaration = "ABC(x1,x2~abc,x3!0,x4~a!1)"
        parsed = parse_molecule(declaration)
        assert parsed["name"] == 'ABC'
        components = parsed["components"]
        assert len(components) == 4

        assert components[0]["name"] == 'x1'
        assert components[0]["state"] == ''
        assert components[0]["bond"] == ''

        assert components[1]["name"] == 'x2'
        assert components[1]["state"] == 'abc'
        assert components[1]["bond"] == ''

        assert components[2]["name"] == 'x3'
        assert components[2]["state"] == ''
        assert components[2]["bond"] == '0'

        assert components[3]["name"] == 'x4'
        assert components[3]["state"] == 'a'
        assert components[3]["bond"] == '1'

    @pytest.mark.parametrize("declaration", [
        "A( )", "ABC( )", "A ()"
        "A(x~a, x~a,y)", "A(x1~a!0, x1~b ,y!1)", "ABC(x1,x2~abc,x3!0 )", "C( x)", "DD(y )"
    ])
    def test_invalid(self, declaration: str):
        with pytest.raises(ValueError):
            parse_molecule(declaration)


class TestParsePattern():
    @pytest.mark.parametrize("declaration, molecules_count", [
        ("A()", 1),
        ("ABC()", 1),
        ("A() ", 1),
        ("A(x~a,x~a,y)", 1),
        ("A(x1~a!0,x1~b,y!1)", 1),
        ("ABC(x1,x2~abc,x3!0,x4~a!1)", 1),
    ])
    def test_single_molecule(self, declaration: str, molecules_count: int):
        parsed = parse_pattern(declaration)
        assert len(parsed) == molecules_count

    @pytest.mark.parametrize("declaration, expected_name, expected_components", [
        ("ABC(x1,x2~abc,x3!0,x4~a!1)", "ABC", [
         ('x1', '', ''), ('x2', 'abc', ''), ('x3', '', '0'), ('x4', 'a', '1')])
    ])
    def test_single_molecule_components(self,
                                        declaration: str,
                                        expected_name: str,
                                        expected_components: list[tuple[str, str, str]]):
        declaration = "ABC(x1,x2~abc,x3!0,x4~a!1)"
        parsed = parse_pattern(declaration)

        assert len(parsed) == 1
        parsed_molecule = parsed[0]
        assert parsed_molecule["name"] == expected_name

        components = parsed_molecule["components"]
        assert len(components) == len(expected_components)

        for comp, expected_comp in zip(components, expected_components):
            assert comp['name'] == expected_comp[0]
            assert comp['state'] == expected_comp[1]
            assert comp['bond'] == expected_comp[2]

    @pytest.mark.parametrize("declaration, molecules_count", [
        ("L(rec!0).R(lig!0)", 2),
        ("R(lig!0).L(rec!0)", 2),
        ("L(rec!1).R(lig!1)", 2),
        ("R(lig!1).L(rec!1)", 2),
        ("L(rec!0).R(lig!0,ch~open)", 2),
        ("L(rec!0).R(lig!0,ch~closed)", 2),
        ("R(lig!0,ch~open).L(rec!0)", 2),
        ("R(lig!0,ch~closed).L(rec!0)", 2),
        ("A(x,y!0).B(p!0,q~a!1).C(r~d!1,s!2).A(x!2,y)", 4),
        ("B(p!0,q~a!1).C(r~d!1,s!2).A(x!2,y).A(x,y!0)", 4),
        ("L(rec!0).R(lig!0) ", 2),
        (" L(rec!0).R(lig!0)", 2)
    ])
    def test_complexes(self, declaration: str, molecules_count: int):
        parsed = parse_pattern(declaration)
        assert len(parsed) == molecules_count

    # @pytest.mark.parametrize("declaration", [
    #     "L( rec!0).R(lig!0)",
    #     "R(lig!0).L(rec!0 )",
    #     "L(rec!1) .R(lig!1)",
    #     "R(lig!1). L(rec!1)",
    #     "L(rec!0 ).R(lig!0,ch~open)",
    #     "L(rec!0).R( lig!0,ch~closed)",
    # ])
    # def test_invalid_complexes(self, declaration: str):
    #     with pytest.raises(ValueError):
    #         parse_pattern(declaration)

    def test_complex_1(self):
        declaration = "A(x,y!0).B(p!0,q~a!1).C(r~d!1,s!2).A(x!2,y)"
        parsed = parse_pattern(declaration)
        assert len(parsed) == 4
        assert parsed[0]["name"] == 'A'
        assert parsed[1]["name"] == 'B'
        assert parsed[2]["name"] == 'C'
        assert parsed[3]["name"] == 'A'
        assert len(parsed[0]["components"]) == 2
        assert len(parsed[1]["components"]) == 2
        assert len(parsed[2]["components"]) == 2
        assert len(parsed[3]["components"]) == 2


class TestParseReaction():
    @pytest.mark.parametrize("declaration, count1, count2", [
        ("R1: S(x~u) + A() + B() + C() -> S(x~p) + A() + B() k1", 4, 3),
        ("R2: S(x~u) + A(b~u!0).B(a!0,c!1).C(b~abc!1) + D(a!0).A(d!0) -> S(x~u) + D(a!0).A(d!0) k2", 3, 2),
        ("R3: R(lig,ch~closed) + L(rec) <-> R(lig!0,ch~closed).L(rec!0) k2, kr", 2, 1)
    ])
    def test_valid_unidirectional(self, declaration: str, count1: int, count2: int):
        parsed = parse_reaction_rule(declaration)
        assert len(parsed["reactants"]) == count1
        assert len(parsed["products"]) == count2


class TestParseReactants():
    @pytest.mark.parametrize("declaration, count", [
        ("S(x~u) + A() + B() + C()", 4),
        ("S(x~u) + A(b~u!0).B(a!0,c!1).C(b~abc!1) + D(a!0).A(d!0)", 3)
    ])
    def test_valid(self, declaration: str, count: int):
        parsed = parse_reactants_sum(declaration)
        assert len(parsed) == count


class TestParseObservable():
    @pytest.mark.parametrize("declaration", [
        "Molecules A A(b)",
        "Species AC A(b,c!1).C(a!0)"
    ])
    def test_valid(self, declaration: str):
        parsed = parse_observable(declaration)
        assert parsed is not None
