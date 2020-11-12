import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.sans-serif": ["Helvetica"],
    "font.family": "DejaVu Sans"})

# Time conversion
second_to_day = 1/(60*60*24)
day_to_second = 60*60*24

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


def setupMatrixCN(n, s):
    # Setup the matrix A
    A = np.zeros([n-1, n-1])
    A[0, 0] = 1 - s/2
    A[0, 1] = s/2
    A[n-2, n-2] = 1 - s/2
    A[n-2, n-3] = s/2

    B = np.zeros([n-1, n-1])
    B[0, 0] = s/2 + 1
    B[0, 1] = -s/2
    B[n-2, n-2] = s/2 + 1
    B[n-2, n-3] = -s/2

    for i in range(1, n-2):
        A[i, i] = 1 - s
        A[i, i - 1] = s/2
        A[i, i + 1] = s/2
        B[i, i] = 1 + s
        B[i, i - 1] = - s/2
        B[i, i + 1] = - s/2
    return np.matmul(np.linalg.inv(B), A)


def solverFTCS(s, n, rho, e, tEnd, Sj, aj, dt):
    A = setupMatrixFTCS(n, s)
    # Forcing constant term
    F = q * Sj * (1-aj)*dt

    t = 0
    T = np.zeros(n-1)
    while t < tEnd:
        t += dt
        T = np.dot(A, T) + F - Q*(T**4) * dt
    return T


def solverCN(s, n, rho, e, tEnd, Sj, aj, dt):
    A = setupMatrixCN(n, s)
    # Forcing constant term
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
n = 182  # Integration points

Pole_albedo = 0.7  # Albedo at the poles as solution to a - b*cos(theta).
Ecvator_albedo = 0.45  # Set to get the albedo at evcator

# Semi fixed variables
startThe = -np.pi/2  # start position in radian
endThe = np.pi/2  # end position in radian
dtheta = np.pi / 180  # theta[1] - theta[0] Steplength in radian
t_end = 140 * 10 * day_to_second  # 10*T_rad in seconds, end time

# Fixed variables
rho = 1  # kg/m^3
heat_cap = 1004  # j/kgK heat cap of air
H = 1e4  # m Height of the atmosphere
sigma = 5.67e-8  # W/m^2k
radiusE = 6.37e6  # m

# Defining fixed constants for the rest of the problem.
q = 1 / (rho * heat_cap * H)
Q = (epsilon * sigma) / (rho * heat_cap * H)

"""Following the the code that performs the calculation and creates the plots"""
if __name__ in "__main__":
    # Setup the problem
    theta = np.linspace(startThe, endThe, n-1)
    Sj, aj = setupArrays(Pole_albedo, Ecvator_albedo, n, theta)

    Dlist = [1e7]  # m^2/s

    for D in Dlist:
        s = 0.4
        # Setup non-constant constants
        dt = s * (radiusE * radiusE * dtheta * dtheta)/D  # time step

        # Solve the problem
        T_ftcs1 = solverFTCS(s, n, rho, epsilon, t_end, Sj, aj, dt)

        T_cn1 = solverCN(s, n, rho, epsilon, t_end, Sj, aj, dt)

        s = 4
        # Setup non-constant constants
        dt = s * (radiusE * radiusE * dtheta * dtheta)/D  # time step

        # Solve the problem
        # T_ftcs2 = solverFTCS(s, n, rho, epsilon, t_end, Sj, aj, dt)

        T_cn2 = solverCN(s, n, rho, epsilon, t_end, Sj, aj, dt)

    plt.plot(np.rad2deg(theta), (T_ftcs1-273), 'k', linewidth=2,
             label=rf"Plot for s = 0.4 with FTCS")

    plt.plot(np.rad2deg(theta), (T_cn1-273), 'r-.',
             label=rf"Plot for s = 0.4 with CN")

    plt.xlabel(r"Longitude $[^{\circ}\theta]$")
    plt.ylabel(r"Temperature $[^{\circ}C]$")
    plt.title(r"Comparison between FTCS and CN")
    plt.legend()
    plt.show()

    plt.plot(np.rad2deg(theta), (T_cn1-273), 'k',
             label=rf"Plot for s = 0.4 with CN")

    plt.plot(np.rad2deg(theta), (T_cn2-273), 'r--',
             label=rf"Plot for s = 4 with CN")

    plt.xlabel(r"Longitude $[^{\circ}\theta]$")
    plt.ylabel(r"Temperature $[^{\circ}C]$")
    plt.title(r"Comparison of CN for different s.")
    plt.legend()
    plt.show()
