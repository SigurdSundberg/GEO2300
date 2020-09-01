import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("earth_density_profile.txt")
z = data[:, 0]  # m
rho = data[:, 1]  # kg/m3


plt.plot(rho, z)
plt.ylabel("[m]")
plt.xlabel("kg/m3")
plt.axvline(x=3.37e3)
plt.axvline(x=5.59e3)
plt.axvline(x=9.86e3)
plt.axvline(x=1.22e4)
plt.axvline(x=1.27e4)
plt.show()
