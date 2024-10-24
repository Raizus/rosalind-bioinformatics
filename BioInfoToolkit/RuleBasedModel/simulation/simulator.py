import abc
import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import write_data_row
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict


class SimulatorABC(abc.ABC):
    sim_params: SimulateDict
    cdat_filename: str = 'output.cdat'
    gdat_filename: str = 'output.gdat'
    next_recording_idx: int

    def __init__(self,
                 sim_params: SimulateDict,
                 cdat_filename: str = 'output.cdat',
                 gdat_filename: str = 'output.gdat'
                 ) -> None:
        super().__init__()
        self.sim_params = sim_params
        self.cdat_filename = cdat_filename
        self.gdat_filename = gdat_filename
        self.next_recording_idx = 0

    def record_line(self,
                    time: float,
                    y: npt.NDArray[np.float_],
                    weights: npt.NDArray[np.float_],
                    recording_times: npt.NDArray[np.float_] | None):
        n_steps = self.sim_params['n_steps']

        # If n_steps is provided, only record concentrations at predefined time points
        if n_steps is not None and recording_times is not None:
            while (self.next_recording_idx < n_steps
                   and time >= recording_times[self.next_recording_idx]):

                groups_conc = y @ weights

                # Record concentrations of reactants
                recording_time = recording_times[self.next_recording_idx]
                write_data_row(self.cdat_filename, recording_time, y)

                # Record concentrations of observables
                write_data_row(self.gdat_filename,
                                recording_time, groups_conc)

                print(f"\tt = {recording_time}")
                self.next_recording_idx += 1
        else:
            # Record time and concentrations for all reactions
            groups_conc = y @ weights

            # Record concentrations of reactants
            write_data_row(self.cdat_filename, time, y)

            # Record concentrations of observables
            write_data_row(self.gdat_filename, time, groups_conc)

            print(f"\tt = {time}")
