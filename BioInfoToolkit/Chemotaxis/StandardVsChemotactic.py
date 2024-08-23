import numpy as np
import matplotlib.pyplot as plt

from BioInfoToolkit.Chemotaxis.CellState import CellState, tumble_move_normal
from BioInfoToolkit.Chemotaxis.Chemotaxis import BacteriumMotionSim, ligand_concentration_2


def interpolate_positions(times_v: np.ndarray, times: np.ndarray, positions: np.ndarray):
    interpolated_positions_x = np.interp(
        times_v, times, positions[:, 0])
    interpolated_positions_y = np.interp(
        times_v, times, positions[:, 1])
    interpolated_positions = np.vstack(
        (interpolated_positions_x, interpolated_positions_y)).T
    return interpolated_positions


def distance_to_goal(n: int):
    total_time = 800
    sim = BacteriumMotionSim(CellState(), total_time=total_time)
    num_time_points = 1000  # Number of points in the time vector
    time_vector = np.linspace(0, total_time, num_time_points)

    # Arrays to store interpolated distances
    all_distances_standard = np.zeros((n, num_time_points))
    all_distances_chemotactic = np.zeros((n, num_time_points))

    # Run n simulations for both strategies
    for i in range(n):
        positions_standard, times_standard = sim.standard_random_walk()
        positions_chemotactic, times_chemotactic = sim.chemotactic_random_walk()

        # Interpolate positions to the common time vector
        interpolated_positions_standard = interpolate_positions(
            time_vector, times_standard, positions_standard)

        interpolated_positions_chemotactic = interpolate_positions(
            time_vector, times_chemotactic, positions_chemotactic)

        # Calculate distances to the goal for each position
        distances_standard = np.linalg.norm(
            interpolated_positions_standard - sim.goal, axis=1)
        distances_chemotactic = np.linalg.norm(
            interpolated_positions_chemotactic - sim.goal, axis=1)

        # Store distances
        all_distances_standard[i, :] = distances_standard
        all_distances_chemotactic[i, :] = distances_chemotactic

    # Calculate mean and standard deviation
    mean_distances_standard = np.mean(all_distances_standard, axis=0)
    std_distances_standard = np.std(all_distances_standard, axis=0)

    mean_distances_chemotactic = np.mean(all_distances_chemotactic, axis=0)
    std_distances_chemotactic = np.std(all_distances_chemotactic, axis=0)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(time_vector, mean_distances_standard,
                color='blue', label='Standard Walk')
    plt.fill_between(time_vector,
                        mean_distances_standard - std_distances_standard,
                        mean_distances_standard + std_distances_standard,
                        color='blue', alpha=0.3)

    plt.plot(time_vector, mean_distances_chemotactic,
                color='green', label='Chemotactic Walk')
    plt.fill_between(time_vector,
                        mean_distances_chemotactic - std_distances_chemotactic,
                        mean_distances_chemotactic + std_distances_chemotactic,
                        color='green', alpha=0.3)

    plt.title('Average Distance to Goal over Time')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Distance to Goal (µm)')
    plt.legend()
    plt.grid(True)
    plt.show()


def mean_run_time_vs_distance(n: int):
    mean_run_times = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
    total_time = 1000
    avg_distances = []
    std_distances_list = []

    num_time_points = 1000  # Number of points in the time vector
    time_vector = np.linspace(0, total_time, num_time_points)

    for mean_run_time in mean_run_times:
        cell_state = CellState(mean_run_time=mean_run_time)
        sim = BacteriumMotionSim(cell_state, total_time=total_time)
        all_distances = np.zeros((n, num_time_points))

        for i in range(n):
            positions, times = sim.chemotactic_random_walk()

            # Interpolate positions to the common time vector
            interpolated_positions = interpolate_positions(time_vector, times, positions)

            # Calculate distances to the goal for each position
            distances = np.linalg.norm(
                interpolated_positions - sim.goal, axis=1)
            all_distances[i, :] = distances

        # Calculate mean and standard deviation of distances
        mean_distances = np.mean(all_distances, axis=0)
        std_distances = np.std(all_distances, axis=0)

        avg_distances.append(mean_distances)
        std_distances_list.append(std_distances)

    # Plotting
    plt.figure(figsize=(12, 8))

    for i, mean_run_time in enumerate(mean_run_times):
        plt.plot(time_vector, avg_distances[i],
                    label=f'µ = {mean_run_time} s')
        plt.fill_between(time_vector,
                         avg_distances[i] - std_distances_list[i],
                         avg_distances[i] + std_distances_list[i],
                         alpha=0.3)

    plt.title('Average Distance to Goal vs. Mean Run Time')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Distance to Goal (µm)')
    plt.legend()
    plt.grid(True)
    plt.show()


def tumble_orientation_vs_distance(n: int):
    total_time = 800
    num_time_points = 1000  # Number of points in the time vector

    cell_state1 = CellState()
    sim1 = BacteriumMotionSim(cell_state1, total_time=total_time)
    cell_state2 = CellState(tumble_move_fun=tumble_move_normal)
    sim2 = BacteriumMotionSim(cell_state2, total_time=total_time,)
    time_vector = np.linspace(0, total_time, num_time_points)

    # Arrays to store interpolated distances
    all_distances_standard = np.zeros((n, num_time_points))
    all_distances_chemotactic = np.zeros((n, num_time_points))

    # Run n simulations for both strategies
    for i in range(n):
        positions_uniform, times_uniform = sim1.chemotactic_random_walk()
        positions_normal, times_normal = sim2.chemotactic_random_walk()

        # Interpolate positions to the common time vector
        interpolated_positions_uniform = interpolate_positions(
            time_vector, times_uniform, positions_uniform)

        interpolated_positions_normal = interpolate_positions(
            time_vector, times_normal, positions_normal)

        # Calculate distances to the goal for each position
        distances_standard = np.linalg.norm(
            interpolated_positions_uniform - sim1.goal, axis=1)
        distances_chemotactic = np.linalg.norm(
            interpolated_positions_normal - sim1.goal, axis=1)

        # Store distances
        all_distances_standard[i, :] = distances_standard
        all_distances_chemotactic[i, :] = distances_chemotactic

    # Calculate mean and standard deviation
    mean_distances_standard = np.mean(all_distances_standard, axis=0)
    std_distances_standard = np.std(all_distances_standard, axis=0)

    mean_distances_chemotactic = np.mean(all_distances_chemotactic, axis=0)
    std_distances_chemotactic = np.std(all_distances_chemotactic, axis=0)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(time_vector, mean_distances_standard,
             color='blue', label='Tumble Uniform Distribution')
    plt.fill_between(time_vector,
                     mean_distances_standard - std_distances_standard,
                     mean_distances_standard + std_distances_standard,
                     color='blue', alpha=0.3)

    plt.plot(time_vector, mean_distances_chemotactic,
             color='green', label='Tumble Normal Distribution')
    plt.fill_between(time_vector,
                     mean_distances_chemotactic - std_distances_chemotactic,
                     mean_distances_chemotactic + std_distances_chemotactic,
                     color='green', alpha=0.3)

    plt.title('Chemotactic Walk - Average Distance to Goal over Time')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Distance to Goal (µm)')
    plt.legend()
    plt.grid(True)
    plt.show()


def single_vs_double_goal(n: int):
    total_time = 800
    num_time_points = 1000  # Number of points in the time vector
    time_vector = np.linspace(0, total_time, num_time_points)

    sim1 = BacteriumMotionSim(CellState(), total_time=total_time)
    sim2 = BacteriumMotionSim(CellState(), total_time=total_time, ligand_concentration=ligand_concentration_2)
    goal1 = np.array([1500, 1500])
    goal2 = np.array([-1500, 1500])

    # Arrays to store interpolated distances
    all_distances_single_goal = np.zeros((n, num_time_points))
    all_distances_double_goal = np.zeros((n, num_time_points))

    # Run n simulations for both strategies
    for i in range(n):
        positions_single, times_single = sim1.chemotactic_random_walk()
        positions_double, times_double = sim2.chemotactic_random_walk()

        # Interpolate positions to the common time vector
        interpolated_positions_single = interpolate_positions(
            time_vector, times_single, positions_single)

        interpolated_positions_double = interpolate_positions(
            time_vector, times_double, positions_double)

        # Calculate distances to the goal for each position
        distances_single = np.linalg.norm(
            interpolated_positions_single - goal1, axis=1)
        
        distances_double_1 = np.linalg.norm(
            interpolated_positions_double - goal1, axis=1)
        distances_double_2 = np.linalg.norm(
            interpolated_positions_double - goal2, axis=1)
        distances_double = np.min( np.column_stack((distances_double_1, distances_double_2)), axis=1)

        # Store distances
        all_distances_single_goal[i, :] = distances_single
        all_distances_double_goal[i, :] = distances_double

    # Calculate mean and standard deviation
    mean_distances_single = np.mean(all_distances_single_goal, axis=0)
    std_distances_single = np.std(all_distances_single_goal, axis=0)

    mean_distances_double = np.mean(all_distances_double_goal, axis=0)
    std_distances_double = np.std(all_distances_double_goal, axis=0)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(time_vector, mean_distances_single,
                color='blue', label='Single Goal')
    plt.fill_between(time_vector,
                        mean_distances_single - std_distances_single,
                        mean_distances_single + std_distances_single,
                        color='blue', alpha=0.3)

    plt.plot(time_vector, mean_distances_double,
                color='green', label='Double Goal')
    plt.fill_between(time_vector,
                        mean_distances_double - std_distances_double,
                        mean_distances_double + std_distances_double,
                        color='green', alpha=0.3)

    plt.title('Average Distance to Nearest Goal over Time')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Distance to Goal (µm)')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    # distance_to_goal(n=200)
    # mean_run_time_vs_distance(300)
    # tumble_orientation_vs_distance(200)
    single_vs_double_goal(400)
