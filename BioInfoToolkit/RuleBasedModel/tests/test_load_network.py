
import pytest

from BioInfoToolkit.RuleBasedModel.model.load_model import load_model
from BioInfoToolkit.RuleBasedModel.network.load_network import load_network
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork
from BioInfoToolkit.RuleBasedModel.utils.plotting import plot_concentrations
from BioInfoToolkit.RuleBasedModel.utils.utls import compose_path, decompose_path

@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/test.net',
    './BioInfoToolkit/RuleBasedModel/assets/test2.net'
])
def test_load_network(fp: str):
    network = load_network(fp)
    assert isinstance(network, ReactionNetwork)


@pytest.mark.parametrize("fp", [
    # './BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.bngl',
    # './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/BLBR.bngl',
])
def test_build_network(fp: str):
    model = load_model(fp)
    network = ReactionNetwork()
    max_stoich: dict[str, int] = {
        'R': 5,
        'L': 5
    }
    network.build_network(model, max_stoich=max_stoich)
    path, name, ext = decompose_path(fp)
    net_path = compose_path(path, f"{name}_gen", ext)
    network.save_network(net_path)
    a = 0


@pytest.mark.parametrize("fp", [
    # './BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.bngl'
])
def test_simulation_network(fp: str):
    model = load_model(fp)
    network = ReactionNetwork()
    network.build_network(model)
    path, name, ext = decompose_path(fp)
    net_path = compose_path(path, f"{name}_gen", ext)
    network.save_network(net_path)
    times, concent = network.gillespie_simulation(1.0, 3000)
    plot_concentrations(times, concent, network.groups_block.items)
