import random
import matplotlib.animation as animation
from time import sleep
import numpy as np
import scipy.signal
# import nptyping as npt
import scipy
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
# npt.NDArray[npt.Shape["Size, Size"], npt.Int]

def diffuse(concentrations: np.ndarray, d: float) -> np.ndarray:
    shape = concentrations.shape
    dtype = concentrations.dtype
    out = np.zeros(shape, dtype=dtype)

    kernel = np.array([
        [0.05*d, 0.2*d, 0.05*d],
        [0.2*d, 1-d, 0.2*d],
        [0.05*d, 0.2*d, 0.05*d]])

    out = scipy.signal.convolve2d(concentrations, kernel, mode='same')

    return out


class GrayScottModel:
    f: float = 0.3     # feed rate
    k: float = 0.4     # kill rate
    dA: float = 0.2    # diffusion rate of A
    dB: float = 0.1    # diffusion rate of B
    r:  float = 1.0    # rate of the reproduction reaction A + 2B -> 3B
    Ac: np.ndarray     # concentrations of A particles
    Bc: np.ndarray     # concentrations of B particles

    def __init__(self,
                 grid_shape: tuple[int, int],
                 f: float = 0.3, 
                 k: float = 0.4,
                 dA: float = 0.2,
                 dB: float = 0.1,
                 r:  float = 1.0) -> None:
        self.f = f
        self.k = k
        self.dA = dA
        self.dB = dB
        self.r = r

        # initial concentrations of A and B particles
        h, w = grid_shape
        Ac = np.ones(grid_shape, dtype=np.float64)
        Bc = np.zeros(grid_shape, dtype=np.float64)

        # Bc[h//2, w//2] = 1

        self.Ac = Ac
        self.Bc = Bc


    def update(self):
        C = self.r * self.Ac * (self.Bc * self.Bc)
        Ac = diffuse(self.Ac, self.dA) + self.f * (1-self.Ac) - C
        Bc = diffuse(self.Bc, self.dB) - self.k * self.Bc + C

        self.Ac = Ac
        self.Bc = Bc


# Function to calculate the pixel value
def calculate_pixel_value(Ac, Bc):
    total = Ac + Bc
    with np.errstate(divide='ignore', invalid='ignore'):
        value = np.where(total != 0, Bc / total, 0.5)
    return value


def gray_scott_animate():
    # Set up the figure, axis, and plot element
    h, w = 200, 200
    grid_shape = (h, w)  # Define the grid size
    dA = 0.5
    dB = 0.25
    f = 0.044
    k = 0.104
    model = GrayScottModel(grid_shape, f=f, k=k, dA=dA, dB=dB)
    dw, dh = 10, 10
    n = 10

    for _ in range(n):
        y = random.randint(0, h-1-dh)
        x = random.randint(0, w-1-dw)
        model.Bc[x:x+dw, y:y+dw] = 1

    # for _ in range(200):
    #     model.update()

    fig, ax = plt.subplots()
    im = ax.imshow(calculate_pixel_value(model.Ac, model.Bc),
                   cmap='Spectral', interpolation='bilinear')
    
    title = ax.set_title('Iteration: 0')

    # Function to update the frame for animation
    def update(frame):
        model.update()
        im.set_array(calculate_pixel_value(model.Ac, model.Bc))

        # Update the title every 10 iterations
        if frame % 10 == 0:
            title.set_text(f'Iteration: {frame}')

        return [im, title]

    # Create the animation
    nframes = 100000
    interval = 5
    ani = animation.FuncAnimation(
        fig, update, frames=nframes, interval=interval, blit=True)

    # Display the animation
    plt.show()

if __name__ == "__main__":
    gray_scott_animate()
