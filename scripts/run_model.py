
from BioInfoToolkit.RuleBasedModel.actions.actions import apply_actions
from BioInfoToolkit.RuleBasedModel.model.load_model import load_bngl
from BioInfoToolkit.RuleBasedModel.utils.plotting import plot_concentrations
from BioInfoToolkit.RuleBasedModel.utils.utls import read_simulation_data


if __name__ == "__main__":
    fp = "./BioInfoToolkit/RuleBasedModel/assets/chemotaxis1.bngl"
    model, actions = load_bngl(fp)
    network = apply_actions(model, actions)
    if network:
        gdat = network.gdat_filename
        times, concentrations = read_simulation_data(gdat)
        plot_concentrations(times, concentrations)
