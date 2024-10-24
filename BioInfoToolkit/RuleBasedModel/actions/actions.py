
import abc

from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict, \
    parse_generate_network, parse_set_concentration, parse_simulate
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import GenerateNetworkDict


class BNGLACtion(abc.ABC):
    pass

class GenerateNetworkAction(BNGLACtion):
    params: GenerateNetworkDict

    def __init__(self, params: GenerateNetworkDict) -> None:
        self.params = params

    @classmethod
    def from_declaration(cls, declaration: str) -> "GenerateNetworkAction":
        parsed = parse_generate_network(declaration)
        action = GenerateNetworkAction(parsed)
        return action


class SimulateAction(BNGLACtion):
    params: SimulateDict

    def __init__(self, params: SimulateDict) -> None:
        self.params = params

    @classmethod
    def from_declaration(cls, declaration: str) -> "SimulateAction":
        parsed = parse_simulate(declaration)
        action = SimulateAction(parsed)
        return action


class SaveConcentrationsAction(BNGLACtion):
    pass


class ResetConcentrationsAction(BNGLACtion):
    pass


class SetConcentrationAction(BNGLACtion):
    pattern: Pattern
    expression: str

    def __init__(self,
                 pattern: Pattern,
                 expression: str) -> None:
        super().__init__()
        self.pattern = pattern
        self.expression = expression

    @classmethod
    def from_declaration(cls, declaration: str) -> "SetConcentrationAction":
        parsed = parse_set_concentration(declaration)
        expression = parsed['expression']
        parsed_pattern = parsed['pattern']
        pattern = Pattern.from_dict(parsed_pattern, None)
        action = SetConcentrationAction(pattern, expression)
        return action


def apply_actions(model: Model, actions: list[BNGLACtion]):
    network: ReactionNetwork | None = None
    for action in actions:
        if isinstance(action, GenerateNetworkAction):
            network = ReactionNetwork()
            network.generate_network(model, action.params)
            network.save_network(overwrite=action.params['overwrite'])

        elif isinstance(action, SimulateAction) and network:
            network.simulate(action.params)
        elif isinstance(action, SaveConcentrationsAction) and network:
            network.save_concentrations()
        else:
            raise NotImplementedError("Action not implemented.")

    return network
