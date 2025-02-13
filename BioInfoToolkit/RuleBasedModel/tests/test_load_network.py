
import pytest

from BioInfoToolkit.RuleBasedModel.actions.actions import GenerateNetworkAction
from BioInfoToolkit.RuleBasedModel.model.load_model import load_bngl
from BioInfoToolkit.RuleBasedModel.network.reaction_network import load_network
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork
from BioInfoToolkit.RuleBasedModel.network.test_utils import compare_observables_dicts, \
    compare_reactions_dicts, compare_species_dicts

FILES_PATH = './BioInfoToolkit/RuleBasedModel/assets/'

bngl_files = [
    f'{FILES_PATH}chemotaxis1.bngl',
    f'{FILES_PATH}chemotaxis2.bngl',
    f'{FILES_PATH}BLBR.bngl',
    f'{FILES_PATH}repressilator.bngl',
    f'{FILES_PATH}egfr_simple.bngl',
    f'{FILES_PATH}FceRI_ji.bngl',
]

net_files = [
    f'{FILES_PATH}chemotaxis1.net',
    f'{FILES_PATH}chemotaxis2.net',
    f'{FILES_PATH}BLBR.net',
    f'{FILES_PATH}repressilator.net',
    f'{FILES_PATH}egfr_simple.net',
    f'{FILES_PATH}FceRI_ji.net',
]

@pytest.mark.parametrize("fp", net_files)
def test_load_network(fp: str):
    network = load_network(fp)
    assert isinstance(network, ReactionNetwork)


@pytest.mark.parametrize("fp", bngl_files)
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

    obs_res = compare_observables_dicts(network1.groups_block.items,
                                        network2.groups_block.items,
                                        sp_bij_map)

    assert obs_res is True
