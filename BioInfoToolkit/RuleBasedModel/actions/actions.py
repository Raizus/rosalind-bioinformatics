
from BioInfoToolkit.RuleBasedModel.model.Model import Model
from BioInfoToolkit.RuleBasedModel.network.reaction_network import ReactionNetwork
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict, parse_generate_network, parse_simulate
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import GenerateNetworkDict
from BioInfoToolkit.RuleBasedModel.utils.utls import compose_path, decompose_path
import abc


class BNGLACtion(abc.ABC):
    pass

class GenerateNetworkAction(BNGLACtion):
    overwrite: bool | None
    text_reaction: bool | None
    max_stoich: dict[str, int] | None
    max_iter: int | None

    def __init__(self, gen_net_dict: GenerateNetworkDict) -> None:
        self.overwrite = gen_net_dict['overwrite']
        self.max_iter = gen_net_dict['max_iter']
        self.max_stoich = gen_net_dict['max_stoich']
        self.text_reaction = gen_net_dict['text_reaction']

    @classmethod
    def from_declaration(cls, declaration: str) -> "GenerateNetworkAction":
        parsed = parse_generate_network(declaration)
        action = GenerateNetworkAction(parsed)
        return action


class SimulateAction(BNGLACtion):
    method: str
    t_start: float
    t_end: float
    n_steps: int | None
    continue_: bool | None

    def __init__(self, simulate_dict: SimulateDict) -> None:
        self.method = simulate_dict['method']
        self.t_start = simulate_dict['t_start']
        self.t_end = simulate_dict['t_end']
        self.n_steps = simulate_dict['n_steps']
        self.continue_ = simulate_dict['continue_']

    @classmethod
    def from_declaration(cls, declaration: str) -> "SimulateAction":
        parsed = parse_simulate(declaration)
        action = SimulateAction(parsed)
        return action
    
    def as_dict(self) -> SimulateDict:
        params: SimulateDict = {
            'method': self.method,
            't_start': self.t_start,
            't_end': self.t_end,
            'n_steps': self.n_steps,
            'continue_': self.continue_
        }
        return params



def apply_actions(model: Model, actions: list[BNGLACtion]):
    network: ReactionNetwork | None = None
    for action in actions:
        if isinstance(action, GenerateNetworkAction):
            network = ReactionNetwork()
            network.build_network(model, action.max_iter, action.max_stoich)
            network.save_network(overwrite=action.overwrite)

        elif isinstance(action, SimulateAction) and network:
            params = action.as_dict()
            network.simulate(params)

    return network
