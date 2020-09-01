import numpy as np

G = 6.674e-11


def dens_func(r, g):
    vol = 4 / 3 * np.pi * r * r * r
    mass = r * r * g / G
    dens = mass / vol
    print(f"Mass: {mass:.3f}KG, Volume: {vol:.3f}, Density: {dens:5.3f}")


dens_func(6371, 9.81)
dens_func(3389.5, 3.711)
