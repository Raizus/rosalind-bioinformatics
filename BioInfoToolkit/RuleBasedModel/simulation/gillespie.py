from functools import reduce
from typing import OrderedDict
import random

import csv
import numpy as np

from BioInfoToolkit.RuleBasedModel.network.group import ObservablesGroup
from BioInfoToolkit.RuleBasedModel.network.reaction import Reaction
from BioInfoToolkit.RuleBasedModel.utils.utls import write_to_csv


def compute_reaction_rates(
    rate_constants: OrderedDict[int, float],
    concentrations: OrderedDict[int, int],
    reactions: OrderedDict[int, Reaction]
):
    rates: OrderedDict[int, float] = OrderedDict()
    # compute the rates for each reaction
    # the rate is equal to the product of the concentrations of the left
    # side reagents and the reaction rate constant
    for r_id, rxn in reactions.items():
        rate_constant = rate_constants[r_id]
        concent_reactants: list[int] = []
        for reactant in rxn.reactants:
            if reactant in concentrations:
                concent_reactants.append(concentrations[reactant])
            else:
                raise KeyError(
                    f"Species with id '{reactant}' not in the concentrations dictionary.")
        rate = rate_constant * \
            reduce(lambda x, y: x*y, concent_reactants, 1.0)
        rates[r_id] = rate

    return rates


class GillespieSimulator:
    total_time: float
    n_steps: int | None
    cdat_filename: str = 'output.cdat'
    gdat_filename: str = 'output.gdat'

    def __init__(self,
                 total_time: float,
                 n_steps: int | None,
                 cdat_filename: str = 'output.cdat',
                 gdat_filename: str = 'output.gdat'
                ) -> None:
        self.total_time = total_time
        self.n_steps = n_steps
        self.cdat_filename = cdat_filename
        self.gdat_filename = gdat_filename

    def create_gdat_file(self,
                         groups_concentration: OrderedDict[int, list[int]],
                         groups: OrderedDict[int, ObservablesGroup],
                         times: list[float]):
        gdat_filename = self.gdat_filename
        with open(gdat_filename, mode='w', newline='', encoding='utf-8') as gdat_file:
            gdat_writer = csv.writer(gdat_file)

            # Write the header: Time and all observable group IDs
            gdat_writer.writerow(['Time'] +  [groups[key].name for key in groups_concentration.keys()])

            # Write data rows
            for i, t in enumerate(times):
                row = [t] + [groups_concentration[g][i]
                            for g in groups_concentration.keys()]
                gdat_writer.writerow(row)

    def simulate(self,
                 reactions: OrderedDict[int, Reaction],
                 rate_constants: OrderedDict[int, float],
                 concentrations: OrderedDict[int, int],
                 groups: OrderedDict[int, ObservablesGroup]):
        n_steps = self.n_steps
        total_time = self.total_time

        # Initialize observable groups concentrations
        groups_concentration: OrderedDict[int, list[int]] = OrderedDict()
        for g_id, group in groups.items():
            concentration = group.compute_concentration(concentrations)
            groups_concentration[g_id] = [int(concentration)]

        # Apply Gillespie algorithm
        time = 0.0
        times: list[float] = [time]

        # If n_steps is provided, generate time points at which to record data
        if n_steps is not None:
            recording_times = np.linspace(0, total_time, n_steps)
            next_recording_idx = 1  # to track the next recording time index

        # Write the header and first row of .cdat file
        header = ['Time'] + list(concentrations.keys())
        row = [time] + list(concentrations.values())
        write_to_csv(self.cdat_filename, 'w', [header, row])
        r_keys = list(reactions.keys())

        while time < total_time:
            rates = compute_reaction_rates(
                rate_constants, concentrations, reactions)
            total_rate = sum(rates.values())
            if total_rate == 0:
                break

            # Time until the next reaction
            tau = np.random.exponential(1 / total_rate)
            time += tau

            # Determine which reaction occurs
            chosen_reaction = random.choices(
                r_keys, weights=list(rates.values()), k=1)[0]

            # Update concentrations based on the chosen reaction
            rxn = reactions[chosen_reaction]
            for reactant in rxn.reactants:
                concentrations[reactant] -= 1
            for prod in rxn.products:
                concentrations[prod] += 1

            # If n_steps is provided, only record concentrations at predefined time points
            if n_steps is not None:
                while next_recording_idx < n_steps and time >= recording_times[next_recording_idx]:
                    times.append(recording_times[next_recording_idx])

                    # Record concentrations of reactants
                    row = [time] + list(concentrations.values())
                    write_to_csv(self.cdat_filename, 'a', [row])

                    # Update observables concentrations
                    for g_id, group in groups.items():
                        concentration = group.compute_concentration(concentrations)
                        groups_concentration[g_id].append(int(concentration))

                    next_recording_idx += 1
            else:
                # Record time and concentrations for all reactions
                times.append(time)

                row = [time] + list(concentrations.values())
                write_to_csv(self.cdat_filename, 'a', [row])

                # Update observables concentrations
                for g_id, group in groups.items():
                    concentration = group.compute_concentration(concentrations)
                    groups_concentration[g_id].append(int(concentration))

        # Write to .gdat file (observable group concentrations)
        self.create_gdat_file(groups_concentration, groups, times)

        return times, groups_concentration
