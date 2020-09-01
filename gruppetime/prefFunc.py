import numpy as np
import matplotlib.pyplot as plt


def pressure(surface_pressure, fluid_density, g, depth):
    if isinstance(depth, (list, tuple, np.ndarray)):
        pressure = surface_pressure + fluid_density * g * depth
        plt.plot(depth, pressure)
        plt.xlabel("depth [m]")
        plt.ylabel("pressure [_]")
        plt.show()
    return surface_pressure + fluid_density * g * depth

