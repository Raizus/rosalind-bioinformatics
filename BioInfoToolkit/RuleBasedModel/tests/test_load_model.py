
import pytest

from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.model.load_model import load_model


@pytest.mark.parametrize("fp", [
    './BioInfoToolkit/RuleBasedModel/assets/test.bngl',
    './BioInfoToolkit/RuleBasedModel/assets/test2.bngl'
])
def test_load_model(fp: str):
    model = load_model(fp)

    assert isinstance(model, Model)
    assert model.validate() is True
