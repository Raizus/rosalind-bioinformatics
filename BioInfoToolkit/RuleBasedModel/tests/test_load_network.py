
import pytest

from BioInfoToolkit.RuleBasedModel.model.load_model import load_model
from BioInfoToolkit.RuleBasedModel.network.load_network import load_network
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork

@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/test.net',
    './BioInfoToolkit/RuleBasedModel/assets/test2.net'
])
def test_load_network(fp: str):
    network = load_network(fp)
    assert isinstance(network, ReactionNetwork)


@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/test.bngl',
    # './BioInfoToolkit/RuleBasedModel/assets/test2.bngl'
])
def test_build_network(fp: str):
    model = load_model(fp)
    network = ReactionNetwork()
    network.build_network(model)
