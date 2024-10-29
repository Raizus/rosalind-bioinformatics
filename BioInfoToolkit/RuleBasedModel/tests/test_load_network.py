
import pytest

from BioInfoToolkit.RuleBasedModel.actions.actions import GenerateNetworkAction
from BioInfoToolkit.RuleBasedModel.model.load_model import load_bngl
from BioInfoToolkit.RuleBasedModel.network.reaction_network import load_network
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork
from BioInfoToolkit.RuleBasedModel.network.test_utils import compare_reactions_dicts, \
    compare_species_dicts


@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.net',
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.net',
    './BioInfoToolkit/RuleBasedModel/assets/BLBR.net',
    './BioInfoToolkit/RuleBasedModel/assets/repressilator.net',
    './BioInfoToolkit/RuleBasedModel/assets/egfr_simple.net',
    './BioInfoToolkit/RuleBasedModel/assets/FceRI_ji.net',
])
def test_load_network(fp: str):
    network = load_network(fp)
    assert isinstance(network, ReactionNetwork)


@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/BLBR.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/repressilator.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/egfr_simple.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/FceRI_ji.bngl',
])
def test_generate_network(fp: str):
    model, actions = load_bngl(fp)
    for action in actions:
        if isinstance(action, GenerateNetworkAction):
            network = ReactionNetwork()
            network.generate_network(model, action.params)
            network.save_network(overwrite=action.params['overwrite'])


@pytest.mark.parametrize("fp1, fp2", [
    ('./BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.net',
     './BioInfoToolkit/RuleBasedModel/assets/chemotaxis1_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.net',
     './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/BLBR.net',
     './BioInfoToolkit/RuleBasedModel/assets/BLBR_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/repressilator.net',
     './BioInfoToolkit/RuleBasedModel/assets/repressilator_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/egfr_simple.net',
     './BioInfoToolkit/RuleBasedModel/assets/egfr_simple_orig.net'),
    ('./BioInfoToolkit/RuleBasedModel/assets/FceRI_ji.net',
     './BioInfoToolkit/RuleBasedModel/assets/FceRI_ji_orig.net'),
])
def test_compare_networks(fp1: str, fp2: str):
    network1 = load_network(fp1)
    network2 = load_network(fp2)

    sp_bij_map = compare_species_dicts(network1.species_block.items,
                                    network2.species_block.items)
    assert sp_bij_map is not None

    r_bij_map = compare_reactions_dicts(network1.reactions_block.items,
                                        network2.reactions_block.items,
                                        sp_bij_map)

    assert r_bij_map is not None
