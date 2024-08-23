import numpy as np
import matplotlib.pyplot as plt

# L = 1 # ligand molecule concentration
# T = 1 # receptor molecule concentration
# LT = 1 # bound complex
# Forward reaction: L + T -> LT
# Reverse reaction: LT -> L + T

# # Steady State condition: k_bind * [L] * [T] = k_dissociate * [LT]
# # Law of conservation of mass [L] and [T] are always constant across the system
# # [L] + [LT] = l_0 <=> [L] = l_0 - [LT]
# # [T] + [LT] = t_0 <=> [T] = t_0 - [LT]
# # substitute in the steady state equation
# # k_bind * (l_0 - [LT]) * (t_0 - [LT]) = k_dissociate * [LT] <=>
# # k_bind [LT]^2 - (k_bind * l_0 + k_bind * t_0 + k_dissociate) [LT] + l_bind * l_0 * t_0 

class LigandReceptorSim:
    k_bind: float = 0.0146  # um^3 / molecule * s^-1
    k_dissociate: float = 35.0  # s^-1
    l_0: int = 10000  # initial number of ligand molecules
    t_0: int = 7000   # initial number of receptor molecules
    total_time: float = 0.2

    def simulate(self):
        # Initial conditions
        L = self.l_0  # initial number of ligand molecules
        T = self.t_0  # initial number of receptor molecules
        LT = 0   # initial number of bound complexes
        time = 0.0
        k_bind = self.k_bind
        k_dissociate = self.k_dissociate

        # Arrays to store the time evolution of the system
        times = [time]
        L_values = [L]
        T_values = [T]
        LT_values = [LT]

        # Gillespie algorithm simulation
        while time < self.total_time:
            # Calculate reaction rates
            r_bind = k_bind * L * T
            r_dissociate = k_dissociate * LT
            total_rate = r_bind + r_dissociate

            if total_rate == 0:
                break

            # Time until the next reaction
            tau = np.random.exponential(1 / total_rate)
            time += tau

            # Determine which reaction occurs
            if np.random.rand() < r_bind / total_rate:
                # L + T -> LT
                L -= 1
                T -= 1
                LT += 1
            else:
                # LT -> L + T
                L += 1
                T += 1
                LT -= 1

            # Store the results
            times.append(time)
            L_values.append(L)
            T_values.append(T)
            LT_values.append(LT)

        return times, L_values, T_values, LT_values
    
    def steady_state_concentrations(self):
        l_0 = self.l_0  # initial number of ligand molecules
        t_0 = self.t_0  # initial number of receptor molecules
        k_bind = self.k_bind
        k_dissociate = self.k_dissociate

        # Calculate steady state concentrations
        a = k_bind
        b = -(k_bind * l_0 + k_bind * t_0 + k_dissociate)
        c = k_bind * l_0 * t_0

        # Solving the quadratic equation: a*[LT]^2 + b*[LT] + c = 0
        discriminant = b**2 - 4 * a * c
        LT_ss_1 = (-b + np.sqrt(discriminant)) / (2 * a)
        LT_ss_2 = (-b - np.sqrt(discriminant)) / (2 * a)

        # Steady state concentrations
        L_ss_1 = l_0 - LT_ss_1
        L_ss_2 = l_0 - LT_ss_2
        T_ss_1 = t_0 - LT_ss_1
        T_ss_2 = t_0 - LT_ss_2

        if all(x >= 0 for x in [LT_ss_1, L_ss_1, T_ss_1]):
            return L_ss_1, T_ss_1, LT_ss_1
        
        assert all(x >= 0 for x in [LT_ss_2, L_ss_2, T_ss_2])
        return L_ss_2, T_ss_2, LT_ss_2


def plot_system():
    sim = LigandReceptorSim()
    times, L_values, T_values, LT_values = sim.simulate()
    L_ss, T_ss, LT_ss = sim.steady_state_concentrations()

    # Plotting
    plt.figure(figsize=(12, 8))
    plt.plot(times, np.array(L_values) / 1e3, label='[L] (Ligand)', color='blue')
    plt.plot(times, np.array(T_values) / 1e3,
            label='[T] (Receptor)', color='green')
    plt.plot(times, np.array(LT_values) / 1e3, label='[LT] (Complex)', color='red')

    # Plot steady-state lines
    plt.axhline(y=L_ss / 1e3, color='blue', linestyle='--',
                label=f'[L]_ss ≈ {round(L_ss)} Molecules/µm^3')
    plt.axhline(y=T_ss / 1e3, color='green', linestyle='--',
                label=f'[T]_ss ≈ {round(T_ss)} Molecules/µm^3')
    plt.axhline(y=LT_ss / 1e3, color='red', linestyle='--',
                label=f'[LT]_ss ≈ {round(LT_ss)} Molecules/µm^3')

    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (molecules x 10^3)')
    plt.title('Gillespie Algorithm: Ligand-Receptor Dynamics')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    plot_system()