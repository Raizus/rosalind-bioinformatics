from typing import Callable
from matplotlib.cm import get_cmap
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from BioInfoToolkit.Chemotaxis.CellState import CellState

ConcentrationFunc = Callable[[np.ndarray, np.ndarray], float]

def ligand_concentration_1(x: np.ndarray, y: np.ndarray) -> float:
    # Ligand concentration function
    goal = np.array([1500, 1500])
    d = np.sqrt((x - goal[0])**2 + (y - goal[1])**2)
    origin = np.array([0, 0])
    D = float(np.linalg.norm(goal - origin))
    return 100 * 10**(6 * (1 - d / D))


def ligand_concentration_2(x: np.ndarray, y: np.ndarray) -> float:
    # Ligand concentration function
    origin = np.array([0, 0])
    goal1 = np.array([1500, 1500])
    goal2 = np.array([-1500, 1500])

    d1 = np.sqrt((x - goal1[0])**2 + (y - goal1[1])**2)
    D1 = float(np.linalg.norm(goal1 - origin))
    d2 = np.sqrt((x - goal2[0])**2 + (y - goal2[1])**2)
    D2 = float(np.linalg.norm(goal2 - origin))
    return 100 * 10**(6 * (1 - d1 / D1)) + 100 * 10**(6 * (1 - d2 / D2))


def draw_ax(ax: plt.Axes, goal: np.ndarray, origin: np.ndarray, ligand_concentration: ConcentrationFunc):
    ax.set_xlim(-500, 2000)
    ax.set_ylim(-500, 2000)

    # Create a grid for the ligand concentration
    grid_size = 300
    x = np.linspace(-500, 2000, grid_size)
    y = np.linspace(-500, 2000, grid_size)
    X, Y = np.meshgrid(x, y)
    Z = ligand_concentration(X, Y)
    
    # Plot the heatmap
    ax.imshow(Z, extent=(-500, 2000, -500, 2000),
                        origin='lower', cmap='Reds', alpha=0.6)
    
    # Plot the origin and goal
    ax.plot(origin[0], origin[1], 'yo', label='Origin')
    ax.plot(goal[0], goal[1], 'gx', label='Goal')
    

class BacteriumMotionSim:
    goal = np.array([1500, 1500])
    total_time = 800  # Total simulation time in seconds
    cell_state: CellState
    ligand_concentration: ConcentrationFunc

    def __init__(self,
                 cell_state: CellState = CellState(),
                 total_time: float = 800,
                 ligand_concentration: ConcentrationFunc = ligand_concentration_1) -> None:
        self.cell_state = cell_state
        origin = cell_state.origin
        self.goal = np.array([1500, 1500])
        self.D = float(np.linalg.norm(self.goal - origin))
        self.total_time = total_time
        self.ligand_concentration = ligand_concentration


    def standard_random_walk(self):
        self.cell_state.reset()

        total_time = self.total_time
        mean_run_time = self.cell_state.mean_run_time

        position = self.cell_state.position
        time_elapsed = 0.0
        positions = [position.copy()]
        times = [time_elapsed]

        while time_elapsed < total_time:
            # Tumble
            tumble_duration, _, _ = self.cell_state.tumble()

            # Run
            run_duration = np.random.exponential(mean_run_time)
            position = self.cell_state.update_position(run_duration)

            # Update time and positions
            time_elapsed += tumble_duration + run_duration
            positions.append(position.copy())
            times.append(time_elapsed)

        return np.array(positions), np.array(times)


    def chemotactic_random_walk(self):
        self.cell_state.reset()

        total_time = self.total_time
        t_response = self.cell_state.t_response
        t_0 = self.cell_state.mean_run_time
        c = 0.000001

        position = self.cell_state.position
        time_elapsed = 0.0
        positions = [position.copy()]
        times = [time_elapsed]
        previous_concentration = self.ligand_concentration(position[0], position[1])

        tumble_duration, _, _ = self.cell_state.tumble()

        while time_elapsed < total_time:
            current_concentration = self.ligand_concentration(position[0], position[1])
            delta_L = (current_concentration -
                        previous_concentration) / previous_concentration
            previous_concentration = current_concentration

            # Calculate the mean run duration µ
            mu = min(max(t_0 * (1 + 10 * delta_L), c), 4 * t_0)

            # Sample run duration
            p = np.random.exponential(mu)
            if p < t_response:
                position = self.cell_state.update_position(p)
                positions.append(position.copy())
                # tumble
                tumble_duration, _, _ = self.cell_state.tumble()
                time_elapsed += p + tumble_duration

            else:
                position = self.cell_state.update_position(t_response)
                positions.append(position.copy())
                time_elapsed += t_response
                # no tumble

            times.append(time_elapsed)

            if time_elapsed >= total_time:
                break

        return np.array(positions), np.array(times)


    def animate(self, mode='both'):
        origin = self.cell_state.origin
        goal = self.goal

        # Run simulations
        if mode in ['standard', 'both']:
            positions_standard, times_standard = self.standard_random_walk()
        if mode in ['chemotactic', 'both']:
            positions_chemotactic, times_chemotactic = self.chemotactic_random_walk()

        # draw axes
        if mode == 'standard':
            # Plot the bacterium's path
            fig, ax_standard = plt.subplots()
            draw_ax(ax_standard, goal, origin, self.ligand_concentration)
            path_standard, = ax_standard.plot(
                [], [], 'b-', lw=1, label='Standard Path')
            bacterium_standard, = ax_standard.plot(
                [], [], 'ro', label='Bacterium')
            ax_standard.set_title("Standard Random Walk")
        elif mode == 'chemotactic':
            # Plot the bacterium's path
            fig, ax_chemotactic = plt.subplots()
            draw_ax(ax_chemotactic, goal, origin, self.ligand_concentration)
            path_chemotactic, = ax_chemotactic.plot(
                [], [], 'g-', lw=1, label='Chemotactic Path')
            bacterium_chemotactic, = ax_chemotactic.plot(
                [], [], 'ro', label='Bacterium')
            ax_chemotactic.set_title("Chemotactic Random Walk")
        else:  # mode == 'both'
            fig, axs = plt.subplots(1, 2, figsize=(12, 6))
            ax_standard, ax_chemotactic = axs
            draw_ax(ax_standard, goal, origin, self.ligand_concentration)
            draw_ax(ax_chemotactic, goal, origin, self.ligand_concentration)

            # Plot the bacterium's path for both strategies
            path_standard, = ax_standard.plot(
                [], [], 'b-', lw=1, label='Standard Path')
            bacterium_standard, = ax_standard.plot(
                [], [], 'ro', label='Bacterium')

            path_chemotactic, = ax_chemotactic.plot(
                [], [], 'g-', lw=1, label='Chemotactic Path')
            bacterium_chemotactic, = ax_chemotactic.plot(
                [], [], 'ro', label='Bacterium')

            ax_standard.set_title("Standard Random Walk")
            ax_chemotactic.set_title("Chemotactic Random Walk")

        def init():
            if mode == 'standard':
                path_standard.set_data([], [])
                bacterium_standard.set_data([], [])
                return path_standard, bacterium_standard
            elif mode == 'chemotactic':
                path_chemotactic.set_data([], [])
                bacterium_chemotactic.set_data([], [])
                return path_chemotactic, bacterium_chemotactic
            else:  # mode == 'both'
                path_standard.set_data([], [])
                bacterium_standard.set_data([], [])
                path_chemotactic.set_data([], [])
                bacterium_chemotactic.set_data([], [])
                return path_standard, bacterium_standard, path_chemotactic, bacterium_chemotactic

        def update(frame):
            if mode == 'standard':
                path_standard.set_data(
                    positions_standard[:frame, 0], positions_standard[:frame, 1])
                bacterium_standard.set_data(
                    [positions_standard[frame, 0]], [positions_standard[frame, 1]])
                return path_standard, bacterium_standard
            elif mode == 'chemotactic':
                path_chemotactic.set_data(
                    positions_chemotactic[:frame, 0], positions_chemotactic[:frame, 1])
                bacterium_chemotactic.set_data(
                    [positions_chemotactic[frame, 0]], [positions_chemotactic[frame, 1]])
                return path_chemotactic, bacterium_chemotactic
            else:  # mode == 'both'
                all_times = np.union1d(times_standard, times_chemotactic)
                index_standard = np.searchsorted(
                    times_standard, all_times[:frame+1])
                index_chemotactic = np.searchsorted(
                    times_chemotactic, all_times[:frame+1])

                # Interpolate positions
                if index_standard[-1] < len(positions_standard):
                    interpolated_position_standard = positions_standard[index_standard[-1]]
                else:
                    interpolated_position_standard = positions_standard[-1]
                if index_chemotactic[-1] < len(positions_chemotactic):
                    interpolated_position_chemotactic = positions_chemotactic[index_chemotactic[-1]]
                else:
                    interpolated_position_chemotactic = positions_chemotactic[-1]

                path_standard.set_data(
                    positions_standard[:index_standard[-1], 0], positions_standard[:index_standard[-1], 1])
                bacterium_standard.set_data(
                    [interpolated_position_standard[0]], [interpolated_position_standard[1]])

                path_chemotactic.set_data(
                    positions_chemotactic[:index_chemotactic[-1], 0], positions_chemotactic[:index_chemotactic[-1], 1])
                bacterium_chemotactic.set_data(
                    [interpolated_position_chemotactic[0]], [interpolated_position_chemotactic[1]])

                return path_standard, bacterium_standard, path_chemotactic, bacterium_chemotactic

        # Create the animation
        if mode == 'standard':
            ani = animation.FuncAnimation(fig, update, frames=len(
                positions_standard), init_func=init, blit=True, interval=20, repeat=False)
        elif mode == 'chemotactic':
            ani = animation.FuncAnimation(fig, update, frames=len(
                positions_chemotactic), init_func=init, blit=True, interval=20, repeat=False)
        else:  # mode == 'both'
            total_frames = len(np.union1d(times_standard, times_chemotactic))
            ani = animation.FuncAnimation(
                fig, update, frames=total_frames, init_func=init, blit=True, interval=20, repeat=False)

        # Show the plot
        plt.legend()
        plt.show()


    def plot(self, n: int = 10, mode: str = 'both'):
        origin = self.cell_state.origin
        goal = self.goal
        cmap = get_cmap('tab10')

        if mode == 'standard':
            _, ax_standard = plt.subplots()
            draw_ax(ax_standard, goal, origin, self.ligand_concentration)
        elif mode == 'chemotactic':
            _, ax_chemotactic = plt.subplots()
            draw_ax(ax_chemotactic, goal, origin, self.ligand_concentration)
        else:
            _, axs = plt.subplots(1, 2, figsize=(12, 6))
            ax_standard, ax_chemotactic = axs
            draw_ax(ax_standard, goal, origin, self.ligand_concentration)
            draw_ax(ax_chemotactic, goal, origin, self.ligand_concentration)

        for i in range(n):
            if mode in ['standard', 'both']:
                positions_standard, _ = self.standard_random_walk()
                ax_standard.plot(
                    positions_standard[:, 0], positions_standard[:, 1], color=cmap(i % 10), lw=1)
                ax_standard.plot(
                    positions_standard[-1, 0], positions_standard[-1, 1], 'og')
            if mode in ['chemotactic', 'both']:
                positions_chemotactic, _ = self.chemotactic_random_walk()
                ax_chemotactic.plot(
                    positions_chemotactic[:, 0], positions_chemotactic[:, 1], color=cmap(i % 10), lw=1)
                ax_chemotactic.plot(
                    positions_chemotactic[-1, 0], positions_chemotactic[-1, 1], 'og')

        if mode == 'standard':
            ax_standard.set_title("Standard Random Walk")
        elif mode == 'chemotactic':
            ax_chemotactic.set_title("Chemotactic Random Walk")
        else:
            ax_standard.set_title("Standard Random Walk")
            ax_chemotactic.set_title("Chemotactic Random Walk")

        plt.legend()
        plt.show()


if __name__ == "__main__":
    sim = BacteriumMotionSim(CellState())
    sim.animate()
    # sim.plot(3, mode='both')
