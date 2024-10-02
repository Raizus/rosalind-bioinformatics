
import pytest

from BioInfoToolkit.RuleBasedModel.network.load_network import load_network
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork

@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/test.net',
    './BioInfoToolkit/RuleBasedModel/assets/test2.net'
])
def test_load_network(fp: str):
    network = load_network(fp)
    assert isinstance(network, ReactionNetwork)
