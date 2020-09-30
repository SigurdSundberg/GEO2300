import numpy as np
import matplotlib.pyplot as plt
from scipy import special


def diffusion(t, x):
    # A * special.erfc(something)
    return -1


# Could need to scale this value
D = 1  # cm^2 s^-1

t_list = [1, 10, 100]
x = np.linspace(0, 100, 200)

for t in t_list:
    plt.plot(x, diffusion(x, t), label=f"Diffusion for t = {t}")
plt.legend()
plt.show()
