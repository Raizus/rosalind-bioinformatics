
import pytest

from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.model.load_model import load_bngl


@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/chemotaxis2.bngl',
])
def test_load_model(fp: str):
    model, actions = load_bngl(fp)

    assert isinstance(model, Model)
    assert model.validate() is True
