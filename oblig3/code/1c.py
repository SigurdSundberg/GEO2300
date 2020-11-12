import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.sans-serif": ["Helvetica"],
    "font.family": "DejaVu Sans"})

# Time conversion
s_to_day = 1/(60*60*24)
day_to_s = 60*60*24

"""Underneath are the function definitions used"""


def S(t):
    return 350 * np.cos(t) + 150  # W/m^2


def albedo(a, b, t):
    return a - b * np.cos(t)


def T_analytical(S, a, e):
    return (S*(1-a)/(e*sigma))**(1/4)


def setupArrays(a, b, n, t):
    # Setup of S_j and a_j
    S_j = np.zeros(n-1)
    a_j = np.zeros(n-1)
    for i in range(n-1):
        S_j[i] = S(t[i])
        a_j[i] = albedo(a, b, t[i])
    return S_j, a_j


def setupMatrixFTCS(n, s):
    # Setup the matrix A
    A = np.zeros([n-1, n-1])
    A[0, 0] = 1 - s
    A[0, 1] = s
    A[n-2, n-2] = 1 - s
    A[n-2, n-3] = s

    for i in range(1, n-2):
        A[i, i] = 1 - 2*s
        A[i, i - 1] = s
        A[i, i + 1] = s
    return A


def solverFTCS(A, n, rho, e, tEnd, Sj, aj, dt):
    # Forcing constant term
    q = 1 / (rho * heat_cap * H)
    Q = (e * sigma) / (rho * heat_cap * H)

    F = q * Sj * (1-aj)*dt

    t = 0
    T = np.zeros(n-1)
    while t < tEnd:
        t += dt
        T = np.dot(A, T) + F - Q*(T**4) * dt
    return T


""" Following is the variable definitions both fixed and non fixed"""
# Variables
epsilon = 0.6  # Emmisivity
s = 0.4
n = 182  # Integration points


# Semi fixed variables
startThe = -np.pi/2  # start position in radian
endThe = np.pi/2  # end position in radian
dtheta = np.pi / 180  # theta[1] - theta[0] Steplength in radian
t_end = 140 * 10 * day_to_s  # 10*T_rad in seconds, end time

# Fixed variables
rho = 1  # kg/m^3
heat_cap = 1004  # j/kgK heat cap of air
H = 1e4  # m Height of the atmosphere
sigma = 5.67e-8  # W/m^2k
radiusE = 6.37e6  # m


"""Following the the code that performs the calculation and creates the plots"""
if __name__ in "__main__":
    # Setup the problem
    theta = np.linspace(startThe, endThe, n-1)
    A = setupMatrixFTCS(n, s)

    # Normal Albedo
    nPole_albedo = 0.7  # Albedo at the poles as solution to a - b*cos(theta).
    nEcvator_albedo = 0.45  # Set to get the albedo at evcator
    nSj, naj = setupArrays(nPole_albedo, nEcvator_albedo, n, theta)

    # Changed albedo
    Pole_albedo = 0.35  # Albedo at the poles as solution to a - b*cos(theta).
    Ecvator_albedo = 0.1  # Set to get the albedo at evcator
    Sj, aj = setupArrays(Pole_albedo, Ecvator_albedo, n, theta)

    D = 1e7  # m^2/s

    # Setup non-constant constants
    dt = s * (radiusE * radiusE * dtheta * dtheta)/D  # time step

    # Base value
    T = solverFTCS(A, n, rho, epsilon, t_end, Sj, aj, dt)

    plt.plot(np.rad2deg(theta), (T-273), 'g',
             label=rf"Plot for changed albedo.")

    nT = solverFTCS(A, n, rho, epsilon, t_end, nSj, naj, dt)
    plt.plot(np.rad2deg(theta), (nT-273), 'k',
             label=rf"Plot for normal albedo.")

    c_nT = nT - 273.15
    c_T = T - 273.15
    uniform = abs(c_nT - c_T)
    mm = np.max(uniform) - np.min(uniform)
    print(f"The difference in maximum and minimum tempratures is {mm:.3f}")

    plt.annotate(
        rf"Changed $\alpha$: {Pole_albedo:.2f} + {Ecvator_albedo:.2f}$\cos(t)$", (-97, 33.5))
    plt.annotate(
        rf"Normal $\alpha$: {nPole_albedo:.2f} + {nEcvator_albedo:.2f}$\cos(t)$", (-97, 31.9))
    plt.title(r"Distrubution of temperature differences.")
    plt.xlabel(r"Longitude $[\theta^{\circ}]$")
    plt.ylabel(r"Temperature $[C^{\circ}]$")
    plt.legend()
    plt.show()
    plt.title(r"Distrubution of temperature differences.")
    plt.xlabel(r"Longitude $[\theta^{\circ}]$")
    plt.ylabel(r"Temperature $[C^{\circ}]$")
    plt.legend()
    plt.plot(uniform)
    plt.show()
