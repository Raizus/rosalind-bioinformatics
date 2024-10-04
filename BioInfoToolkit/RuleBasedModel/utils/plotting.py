
from typing import OrderedDict
import matplotlib.pyplot as plt
from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup

def plot_concentrations(times: list[float],
                        groups_concentration: OrderedDict[int, list[int]],
                        observables_groups: OrderedDict[int, ObservablesGroup]):
    """Creates a plot of the evolution of the concentrations of different reactants along time

    Args:
        times (list[float]): list of time points
        groups_concentration (OrderedDict[int, list[int]]): dictionary where each element 
        corresponds to a group of reactants and their respective concentration across time
        observables_groups (OrderedDict[int, ObservablesGroup]): dictionary mapping observables 
        index to the object representing the group. Needed for labeling the plot
    """

    plt.figure(figsize=(10, 6))  # Create a figure with a specific size

    for i, concentrations in groups_concentration.items():
        label_i = observables_groups[i].name  # Get the name for labeling
        # Plot each concentration over time
        plt.plot(times, concentrations, label=label_i)

    # Adding title, labels, and legend
    plt.title('Concentration of Reactants Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration')
    # Automatically place the legend at the best location
    plt.legend(loc='best')
    plt.grid(True)  # Enable grid for better readability
    plt.tight_layout()  # Adjust layout to fit everything neatly

    # Show the plot
    plt.show()
    plt.savefig('plot.png')
