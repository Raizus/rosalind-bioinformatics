import abc
import numpy as np
import numpy.typing as npt

from BioInfoToolkit.RuleBasedModel.simulation.simulation_utils import write_data_row
from BioInfoToolkit.RuleBasedModel.utils.action_parsers import SimulateDict
from BioInfoToolkit.RuleBasedModel.utils.utls import get_cdat_last_line


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

    def get_initial(self, y: npt.NDArray[np.float_]) -> tuple[float, npt.NDArray[np.float_]]:
        t_start = self.sim_params['t_start']
        continue_ = self.sim_params['continue_']

        if continue_:
            # read last line
            t, y = get_cdat_last_line(self.cdat_filename)
            if t_start is None or t_start == t:
                t_start = t
            else:
                msg = ("When continue is True, start_time must be "
                       "equal to the last value in the cdat file."
                       f" t_start = {t_start}; t = {t}")
                raise ValueError(msg)

        if t_start is None:
            t_start = 0.0

        return t_start, y

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
