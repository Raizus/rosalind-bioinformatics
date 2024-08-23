import numpy as np
import matplotlib.pyplot as plt

class NegativeAutoregulationMotif:
    k1: float = 1.0      # Rate for X -> X + Y
    k2: float = 0.1      # Rate for Y degradation (Y -> ∅)
    k3: float = 0.01     # Rate for negative autoregulation (2Y -> Y)
    diff_rate: float = 0.05  # Diffusion rate between the cells

    # Simulation parameters
    k1: float = 1.0      # Rate for X -> X + Y
    k2: float = 0.1      # Rate for Y degradation (Y -> ∅)
    k3: float = 0.01     # Rate for negative autoregulation (2Y -> Y)
    diff_rate: float = 0.05  # Diffusion rate between the cells

    # Initial conditions
    X1: int = 100  # Number of X molecules in Cell 1 (steady state)
    X2: int = 100  # Number of X molecules in Cell 2 (steady state)
    Y1: int = 0    # Initial number of Y molecules in Cell 1
    Y2: int = 0    # Initial number of Y molecules in Cell 2
    t: float = 0     # Initial time

    # Time for simulation
    max_time = 50
    time_points = np.linspace(0, max_time, 500)
    Y1_concentration = []
    Y2_concentration = []

    # Function to perform Gillespie step
    def gillespie_step(self, X1: int, X2: int, Y1: int, Y2: int, t: float):
        k1 = self.k1
        k2 = self.k2
        k3 = self.k3
        diff_rate = self.diff_rate
        max_time = self.max_time

        rates = [
            k1 * X1,                # X1 -> X1 + Y1
            k1 * X2,                # X2 -> X2 + Y2
            k2 * Y1,                # Y1 -> ∅ (degradation in Cell 1)
            k2 * Y2,                # Y2 -> ∅ (degradation in Cell 2)
            # 2Y2 -> Y2 (negative autoregulation in Cell 2)
            k3 * Y2 * (Y2 - 1) / 2,
            diff_rate * Y1,         # Y1 diffuses to Y2
            diff_rate * Y2          # Y2 diffuses to Y1
        ]
        total_rate = sum(rates)

        if total_rate == 0:
            return t + max_time, X1, X2, Y1, Y2

        # Time step
        dt = np.random.exponential(1 / total_rate)

        # Which reaction occurs?
        r = np.random.uniform(0, total_rate)
        if r < rates[0]:
            Y1 += 1
        elif r < rates[0] + rates[1]:
            Y2 += 1
        elif r < rates[0] + rates[1] + rates[2]:
            Y1 -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3]:
            Y2 -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4]:
            Y2 -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5]:
            Y1 -= 1
            Y2 += 1
        else:
            Y2 -= 1
            Y1 += 1

        return t + dt, X1, X2, Y1, Y2

    def simulate(self):
        t = self.t
        X1 = self.X1
        X2 = self.X2
        Y1 = self.Y1
        Y2 = self.Y2

        max_time = self.max_time
        Y1_concentration = self.Y1_concentration
        Y2_concentration = self.Y2_concentration

        # Simulation loop
        while t < max_time:
            t, X1, X2, Y1, Y2 = self.gillespie_step(X1, X2, Y1, Y2, t)

            Y1_concentration.append(Y1)
            Y2_concentration.append(Y2)

    def plot_results(self):
        max_time = self.max_time
        Y1_concentration = self.Y1_concentration
        Y2_concentration = self.Y2_concentration
    
        # Plot the results
        plt.figure(figsize=(10, 5))
        plt.plot(np.linspace(0, max_time, len(Y1_concentration)),
                Y1_concentration, label="Cell 1 (no autoregulation)")
        plt.plot(np.linspace(0, max_time, len(Y2_concentration)),
                Y2_concentration, label="Cell 2 (with autoregulation)")
        plt.xlabel("Time")
        plt.ylabel("Concentration of Y")
        plt.legend()
        plt.title("Reaction-Diffusion Model of Transcriptional Regulation")
        plt.grid(visible=True)
        plt.show()


class IncoherentFeedforwardLoopMotif:
    k1: float = 1.0 # Rate for X -> X + Z in Cell 1
    k2: float = 1.5 # Rate for X -> X + Z in Cell 2 (higher rate due to FFL)
    k3: float = 1.0 # Rate for X -> X + Y in Cell 2
    k4: float = 0.001 # Rate for Y + Z -> Y in Cell 2 (repression)
    k_degradation: float = 0.1  # Degradation rate for Y and Z
    diff_rate: float = 0.05  # Diffusion rate between the cells

    # Initial conditions
    X1: int = 100  # Number of X molecules in Cell 1 (steady state)
    X2: int = 100  # Number of X molecules in Cell 2 (steady state)
    Y2: int = 0    # Initial number of Y molecules in Cell 2
    Z1: int = 0    # Initial number of Z molecules in Cell 1
    Z2: int = 0    # Initial number of Z molecules in Cell 2
    t: float = 0   # Initial time

    # Time for simulation
    max_time = 50
    time_points = np.linspace(0, max_time, 500)
    Z1_concentration = []
    Z2_concentration = []

    # Function to perform Gillespie step
    def gillespie_step(self, X1: int, X2: int, Y2: int, Z1: int, Z2: int, t: float):
        k1 = self.k1
        k2 = self.k2
        k3 = self.k3
        k4 = self.k4
        k_degradation = self.k_degradation
        diff_rate = self.diff_rate
        max_time = self.max_time

        rates = [
            k1 * X1,                # X1 -> X1 + Z1 (Cell 1)
            k2 * X2,                # X2 -> X2 + Z2 (Cell 2)
            k3 * X2,                # X2 -> X2 + Y2 (Cell 2)
            k4 * Y2 * Z2,           # Y2 + Z2 -> Y2 (repression in Cell 2)
            k_degradation * Z1,     # Z1 -> ∅ (degradation in Cell 1)
            k_degradation * Z2,     # Z2 -> ∅ (degradation in Cell 2)
            k_degradation * Y2,     # Y2 -> ∅ (degradation in Cell 2)
            diff_rate * Z1,         # Z1 diffuses to Z2
            diff_rate * Z2          # Z2 diffuses to Z1
        ]
        total_rate = sum(rates)

        if total_rate == 0:
            return t + max_time, X1, X2, Y2, Z1, Z2

        # Time step
        dt = np.random.exponential(1 / total_rate)

        # Which reaction occurs?
        r = np.random.uniform(0, total_rate)
        if r < rates[0]:
            Z1 += 1
        elif r < rates[0] + rates[1]:
            Z2 += 1
        elif r < rates[0] + rates[1] + rates[2]:
            Y2 += 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3]:
            Z2 -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4]:
            Z1 -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5]:
            Z2 -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5] + rates[6]:
            Y2 -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5] + rates[6] + rates[7]:
            Z1 -= 1
            Z2 += 1
        else:
            Z2 -= 1
            Z1 += 1

        return t + dt, X1, X2, Y2, Z1, Z2

    def simulate(self):
        t = self.t
        X1 = self.X1
        X2 = self.X2
        Y2 = self.Y2
        Z1 = self.Z1
        Z2 = self.Z2

        max_time = self.max_time
        Z1_concentration = self.Z1_concentration
        Z2_concentration = self.Z2_concentration

        # Simulation loop
        while t < max_time:
            t, X1, X2, Y2, Z1, Z2 = self.gillespie_step(X1, X2, Y2, Z1, Z2, t)

            Z1_concentration.append(Z1)
            Z2_concentration.append(Z2)

    def plot_results(self):
        max_time = self.max_time
        Z1_concentration = self.Z1_concentration
        Z2_concentration = self.Z2_concentration

        # Plot the results
        plt.figure(figsize=(10, 5))
        plt.plot(np.linspace(0, max_time, len(Z1_concentration)),
                 Z1_concentration, label="Cell 1 (simple activation)")
        plt.plot(np.linspace(0, max_time, len(Z2_concentration)),
                 Z2_concentration, label="Cell 2 (with FFL)")
        plt.xlabel("Time")
        plt.ylabel("Concentration of Z")
        plt.legend()
        plt.title("Reaction-Diffusion Model of Incoherent Feedforward Loop")
        plt.grid(visible=True)
        plt.show()


class RepressilatorMotif:
    def __init__(self):
        # Reaction rates
        self.k_prod = 1.0        # Rate for I -> I + X, I -> I + Y, I -> I + Z
        self.k_repress = 0.1     # Rate for X + Y -> Y, Y + Z -> Z, Z + X -> X
        self.k_degradation = 0.1  # Degradation rate for X, Y, Z
        self.diff_rate = 0.05    # Diffusion rate for X, Y, Z

        # Initial conditions
        self.I = 50  # Number of I molecules (steady state activator)
        self.X = 10   # Initial number of X molecules
        self.Y = 0    # Initial number of Y molecules
        self.Z = 0    # Initial number of Z molecules
        self.t = 0    # Initial time

        # Time for simulation
        self.max_time = 20
        self.time_points = np.linspace(0, self.max_time, 1000)
        self.X_concentration = []
        self.Y_concentration = []
        self.Z_concentration = []

    def gillespie_step(self, I, X, Y, Z, t):
        rates = [
            self.k_prod * I,            # I -> I + X
            self.k_prod * I,            # I -> I + Y
            self.k_prod * I,            # I -> I + Z
            self.k_repress * X * Y,     # X + Y -> Y (repression)
            self.k_repress * Y * Z,     # Y + Z -> Z (repression)
            self.k_repress * Z * X,     # Z + X -> X (repression)
            self.k_degradation * X,     # X -> ∅ (degradation)
            self.k_degradation * Y,     # Y -> ∅ (degradation)
            self.k_degradation * Z      # Z -> ∅ (degradation)
        ]
        total_rate = sum(rates)

        if total_rate == 0:
            return t + self.max_time, I, X, Y, Z

        # Time step
        dt = np.random.exponential(1 / total_rate)

        # Which reaction occurs?
        r = np.random.uniform(0, total_rate)
        if r < rates[0]:
            X += 1
        elif r < rates[0] + rates[1]:
            Y += 1
        elif r < rates[0] + rates[1] + rates[2]:
            Z += 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3]:
            Y -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4]:
            Z -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5]:
            X -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5] + rates[6]:
            X -= 1
        elif r < rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5] + rates[6] + rates[7]:
            Y -= 1
        else:
            Z -= 1

        return t + dt, I, X, Y, Z

    def simulate(self):
        t = self.t
        I = self.I
        X = self.X
        Y = self.Y
        Z = self.Z

        while t < self.max_time:
            t, I, X, Y, Z = self.gillespie_step(I, X, Y, Z, t)

            self.X_concentration.append(X)
            self.Y_concentration.append(Y)
            self.Z_concentration.append(Z)

    def plot_results(self):
        # # Plot the results

        # Time points for plotting
        time_points = np.linspace(0, self.max_time, len(self.X_concentration))
        
        # Create a figure with 3 subplots (3 rows, 1 column)
        fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

        # Plot X concentration
        axs[0].plot(time_points, self.X_concentration, color='blue')
        axs[0].set_ylabel("X concentration")
        axs[0].set_title("X Concentration Over Time")
        axs[0].grid(visible=True)

        # Plot Y concentration
        axs[1].plot(time_points, self.Y_concentration, color='green')
        axs[1].set_ylabel("Y concentration")
        axs[1].set_title("Y Concentration Over Time")
        axs[1].grid(visible=True)

        # Plot Z concentration
        axs[2].plot(time_points, self.Z_concentration, color='red')
        axs[2].set_ylabel("Z concentration")
        axs[2].set_xlabel("Time")
        axs[2].set_title("Z Concentration Over Time")
        axs[2].grid(visible=True)

        # Adjust layout
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    # model = NegativeAutoregulationMotif()
    model = RepressilatorMotif()
    model.simulate()
    model.plot_results()