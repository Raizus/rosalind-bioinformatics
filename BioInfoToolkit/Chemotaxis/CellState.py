from typing import Callable
import numpy as np


def tumble_move_uniform(mean_tumble_time: float, prev_angle: float):
    tumble_duration = np.random.exponential(mean_tumble_time)
    angle = np.random.uniform(0, 2 * np.pi)
    direction = np.array([np.cos(angle), np.sin(angle)])
    return tumble_duration, angle, direction


def tumble_move_normal(mean_tumble_time: float, prev_angle: float):
    tumble_duration = np.random.exponential(mean_tumble_time)
    mean = 68 * np.pi / 180 + prev_angle
    std = 36 * np.pi / 180
    angle = np.random.normal(mean, std)
    direction = np.array([np.cos(angle), np.sin(angle)])
    return tumble_duration, angle, direction


class CellState:
    t_response: float = 0.5  # seconds
    mean_tumble_time = 0.1  # seconds
    mean_run_time = 1.0     # seconds (mean tumble period)
    speed = 20.0            # µm/s

    origin = np.array([0, 0])
    angle = np.random.uniform(0, 2 * np.pi)
    direction = np.array([np.cos(0), np.sin(0)])
    tumble_move_fun: Callable[[float, float], tuple[float, float, np.ndarray]]

    def __init__(self,
                 t_response: float = 0.5,
                 mean_run_time: float = 1.0,
                 tumble_move_fun: Callable[[float, float], tuple[float, float, np.ndarray]] = tumble_move_uniform) -> None:
        self.t_response = t_response
        self.mean_run_time = mean_run_time
        self.tumble_move_fun = tumble_move_fun
        self.position = self.origin
        self.tumble()

    def tumble(self):
        tumble_duration, angle, direction = self.tumble_move_fun(self.mean_tumble_time, self.angle)
        self.angle = angle
        self.direction = direction
        return tumble_duration, angle, direction
    
    def reset(self):
        self.position = self.origin
        self.angle = np.random.uniform(0, 2 * np.pi)
        self.tumble()
    
    def update_position(self, t: float):
        self.position = self.position + self.direction * self.speed * t
        return self.position

