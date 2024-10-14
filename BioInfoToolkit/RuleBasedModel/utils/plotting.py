
import matplotlib.pyplot as plt


def plot_concentrations(times: list[float], concentrations: dict[str, list[float]]):
    """
    Creates a plot of the evolution of the concentrations of different reactants along time.

    Args:
        times (list[float]): list of time points.
        concentrations (dict[str, list[float]]): dictionary where the keys are the labels for each reactant and the values are lists of concentrations over time.
    """

    plt.figure(figsize=(10, 6))  # Create a figure with a specific size

    for label, conc_values in concentrations.items():
        # Plot each concentration over time with label as the legend
        plt.plot(times, conc_values, label=label)

    # Adding title, labels, and legend
    plt.title('Concentration of Reactants Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration')
    # Automatically place the legend at the best location
    plt.legend(loc='best')
    plt.grid(True)  # Enable grid for better readability
    plt.tight_layout()  # Adjust layout to fit everything neatly

    # Show the plot
    plt.savefig('plot.png')
    plt.show()
