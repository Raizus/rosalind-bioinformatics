
import pytest

from BioInfoToolkit.RuleBasedModel.model.load_model import load_bngl
from BioInfoToolkit.RuleBasedModel.network.load_network import load_network
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork
from BioInfoToolkit.RuleBasedModel.network.test_utils import compare_reactions_dicts, compare_species_dicts
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
    model, actions = load_bngl(fp)
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
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.bngl',
    # './BioInfoToolkit/RuleBasedModel/assets/test2.bngl',
    # './BioInfoToolkit/RuleBasedModel/assets/BLBR.bngl',
])
def test_simulation_network(fp: str):
    model, actions = load_bngl(fp)
    network = ReactionNetwork()
    # max_stoich: dict[str, int] = {
    #     'R': 5,
    #     'L': 5
    # }
    network.build_network(model)
    network.save_network()


@pytest.mark.parametrize("fp1, fp2", [
    ('./BioInfoToolkit/RuleBasedModel/assets/test2.net',
     './BioInfoToolkit/RuleBasedModel/assets/test2_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.net',
     './BioInfoToolkit/RuleBasedModel/assets/chemotaxis1_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.net',
     './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/BLBR.net',
     './BioInfoToolkit/RuleBasedModel/assets/BLBR_orig.net'),
])
def test_compare_networks(fp1: str, fp2: str):
    network1 = load_network(fp1)
    network2 = load_network(fp2)

    bij_map = compare_species_dicts(network1.species_block.items,
                                    network2.species_block.items)
    assert bij_map is not None

    r_bij_map = compare_reactions_dicts(network1.reactions_block.items,
                                        network2.reactions_block.items,
                                        bij_map)

    assert r_bij_map is not None
