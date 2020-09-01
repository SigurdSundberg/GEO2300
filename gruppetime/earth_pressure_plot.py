import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("earth_density_profile.txt")
z = data[:, 0]  # m
rho = data[:, 1]  # kg/m^3


def pressure(z, rho):
    g = 9.81  # m/s²
    n = len(rho)
    p = np.zeros(n)
    for i in range(n-1):
        del_z = z[i+1] - z[i]
        p[i+1] = p[i]-g*rho[i]*del_z
    return p


p = pressure(z, rho)/1e9
plt.plot(p, z)
plt.ylabel("[m]")
plt.xlabel("GPa")
plt.axhline(y=-6.9e4, color='r')
plt.axhline(y=-2.884e6, color='r')
plt.axhline(y=-5.188e6, color='r')
plt.show()
