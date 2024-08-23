import numpy as np
import matplotlib.pyplot as plt

class ChemotaxisSim:
    NaV: float = 6.02e8  # Unit conversion factor (Molecules/µm^3)

    # Molecule counts (initial)
    L0: int = 5000
    T0: int = 7000
    CheY0: int = 20000
    CheZ0: int = 6000

    # Initial concentrations
    L: int = L0
    T_U: int = int(T0 * 0.8)
    T_P: int = int(T0 * 0.2)
    CheY_U: int = int(CheY0 * 0.5)
    CheY_P: int = int(CheY0 * 0.5)
    CheZ: int = CheZ0

    # Reaction rates
    k_lr_bind: float = 8.8e6 / NaV   # Ligand-receptor binding
    k_lr_dis: float = 35.0           # Ligand-receptor dissociation
    k_T_phos: float = 15.0           # Receptor complex autophosphorylation
    k_Y_phos: float = 3.8e6 / NaV    # Receptor complex phosphorylates CheY
    k_Y_dephos: float = 8.6e5 / NaV  # CheZ dephosphorylates CheY

    total_time: float = 3.0  # Simulation time

    def __init__(self, L0: int = 5000) -> None:
        self.L0 = L0

    def simulate(self):
        # Initial conditions
        L = self.L
        T_U = self.T_U
        T_P = self.T_P
        LT_U = 0
        LT_P = 0
        CheY_U = self.CheY_U
        CheY_P = self.CheY_P
        CheZ = self.CheZ
        time = 0.0

        # Arrays to store the time evolution of the system
        times = [time]
        LT_values = [LT_U + LT_P]
        T_P_values = [T_P]
        CheY_P_values = [CheY_P]

        # Gillespie algorithm simulation
        while time < self.total_time:
            # Calculate reaction rates
            r_bind_U = self.k_lr_bind * L * T_U
            r_dis_U = self.k_lr_dis * LT_U
            r_bind_P = self.k_lr_bind * L * T_P
            r_dis_P = self.k_lr_dis * LT_P

            r_T_phos_U = self.k_T_phos * T_U
            r_T_phos_P = self.k_T_phos * 0.2 * LT_U

            r_Y_phos = self.k_Y_phos * T_P * CheY_U
            r_Y_dephos = self.k_Y_dephos * CheZ * CheY_P

            total_rate = (r_bind_U + r_dis_U + r_bind_P + r_dis_P +
                          r_T_phos_U + r_T_phos_P + r_Y_phos + r_Y_dephos)

            if total_rate == 0:
                break

            # Time until the next reaction
            tau = np.random.exponential(1 / total_rate)
            time += tau

            # Determine which reaction occurs
            reaction_choice = np.random.rand() * total_rate
            if reaction_choice < r_bind_U:
                # L + T_U -> LT_U
                L -= 1
                T_U -= 1
                LT_U += 1
            elif reaction_choice < r_bind_U + r_dis_U:
                # LT_U -> L + T_U
                L += 1
                T_U += 1
                LT_U -= 1
            elif reaction_choice < r_bind_U + r_dis_U + r_bind_P:
                # L + T_P -> LT_P
                L -= 1
                T_P -= 1
                LT_P += 1
            elif reaction_choice < r_bind_U + r_dis_U + r_bind_P + r_dis_P:
                # LT_P -> L + T_P
                L += 1
                T_P += 1
                LT_P -= 1
            elif reaction_choice < r_bind_U + r_dis_U + r_bind_P + r_dis_P + r_T_phos_U:
                # T_U -> T_P
                T_U -= 1
                T_P += 1
            elif reaction_choice < (r_bind_U + r_dis_U + r_bind_P + r_dis_P +
                                    r_T_phos_U + r_T_phos_P):
                # LT_U -> LT_P
                LT_U -= 1
                LT_P += 1
            elif reaction_choice < (r_bind_U + r_dis_U + r_bind_P + r_dis_P +
                                    r_T_phos_U + r_T_phos_P + r_Y_phos):
                # T_P + CheY_U -> T_U + CheY_P
                T_P -= 1
                T_U += 1
                CheY_U -= 1
                CheY_P += 1
            else:
                # CheZ + CheY_P -> CheZ + CheY_U
                CheY_P -= 1
                CheY_U += 1

            # Store the results
            times.append(time)
            LT_values.append(LT_U + LT_P)
            T_P_values.append(T_P)
            CheY_P_values.append(CheY_P)

        return times, LT_values, T_P_values, CheY_P_values


def plot_system():
    sim = ChemotaxisSim()
    times, LT_values, T_P_values, CheY_P_values = sim.simulate()

    # Plotting
    plt.figure(figsize=(12, 8))
    plt.plot(times, np.array(LT_values) / 1e3,
             label='[LT] (Bound Ligand)', color='blue')
    plt.plot(times, np.array(T_P_values) / 1e3,
             label='[T_P] (Phosphorylated CheA)', color='green')
    plt.plot(times, np.array(CheY_P_values) / 1e3,
             label='[CheY_P] (Phosphorylated CheY)', color='red')

    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (molecules x 10^3)')
    plt.title('Gillespie Algorithm: Chemotaxis Signaling Dynamics')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    plot_system()
