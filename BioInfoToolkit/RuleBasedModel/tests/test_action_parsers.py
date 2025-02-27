"""
    Tests for parsing action strings
"""

import pytest
from BioInfoToolkit.RuleBasedModel.utils.action_parsers \
    import parse_generate_network, parse_set_concentration, parse_simulate


class TestParseGenerateNetwork:
    """
        Tests parsing of generate_network actions
    """
    @pytest.mark.parametrize("declaration", [
        "generate_network({overwrite=>1,max_stoich=>{R=>5,L=>5},max_iter=>20,TextReaction=>1});",
        "generate_network({max_iter=>20,TextReaction=>1,overwrite=>1,max_stoich=>{R=>5,L=>5}})",
        "generate_network({TextReaction=>1,max_stoich=>{R=>5}});",
        "generate_network();",
    ])
    def test_valid(self, declaration: str):
        """Test valid generate_network actions

        Args:
            declaration (str): _description_
        """
        parsed = parse_generate_network(declaration)
        assert parsed is not None


class TestParseSimulate:
    """
        Tests parsing of simulate actions
    """
    @pytest.mark.parametrize("declaration, method, t_start, t_end, n_steps", [
        ('simulate({method=>"ssa", t_end=>1, n_steps=>100});',
         ("ssa", None, 1.0, 100)),
        ('simulate({method=>"ssa", t_end=>0.8})', ("ssa", None, 0.8, None)),
        ('simulate({method=>"ssa", t_start=>0.8, t_end=>4.5})',
         ("ssa", 0.8, 4.5, None)),
    ])
    def test_valid(self, declaration: str,
                   params: tuple[str, float | None, float | None, int | None]):
        """Test valid simulate actions

        Args:
            declaration (str): _description_
            params (tuple[str, float|None, float|None, int|None]): action parameters
        """
        parsed = parse_simulate(declaration)
        assert params[0] == parsed["method"]
        assert params[1] == parsed["t_start"]
        assert params[2] == parsed["t_end"]
        assert params[3] == parsed["n_steps"]


class TestParseSetConcentration:
    """
        Tests parsing of setConcentration actions
    """
    @pytest.mark.parametrize("declaration, expression", [
        ('setConcentration("Lig(l,l)", "Lig_tot")', "Lig_tot"),
    ])
    def test_valid(self, declaration: str, expression: str):
        parsed = parse_set_concentration(declaration)
        assert parsed["expression"] == expression
