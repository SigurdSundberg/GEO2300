import matplotlib.pyplot as plt
import numpy as np
from math import ceil

plt.rcParams.update({
    "text.usetex": True,
    "font.sans-serif": ["Helvetica"],
    "font.family": "DejaVu Sans"})


def pdf(x, sigma, mu):
    leading_term = 1/(sigma * np.sqrt(2*np.pi))
    divisor = (2*sigma*sigma)
    return leading_term * np.exp(-(x-mu)*(x-mu) / divisor)


def analytical(x, t, u, sigma, mu):
    leading_term = 1/(sigma * np.sqrt(2*np.pi))
    divisor = (2*sigma*sigma)
    return leading_term * np.exp(-(x-u*t-mu)*(x-u*t-mu) / divisor)


def FTCS_matrix(n, c):
    A = np.zeros([n, n])
    upper = - c / 2
    lower = c / 2
    for i in range(1, n-1):
        A[i, i] = 1
        A[i, i-1] = lower
        A[i, i+1] = upper
    A[0, 0] = 1
    A[-1, -1] = 1
    A[0, 1] = upper
    A[-1, -2] = lower
    return A


if __name__ == "__main__":
    end_times = [0, 500, 1000, 1500]  # [s]
    C = 0.1  # courant number
    mu = 100  # [m] Relative mean
    sigma = 10  # [m] Standard deviation

    dx = 5
    n = ceil(1000/dx)
    x = np.linspace(0, 1000, n+1)  # meters
    u = 0.5  # m/s

    dt = C/u * dx

    A = FTCS_matrix(n+1, C)

    stability = C * dt / dx
    print(f"Stability r: {stability:.3f}")

    for t_end in end_times:
        t = 0
        p = pdf(x, sigma, mu)
        while t < t_end:
            t += dt
            p = np.dot(A, p)
        plt.plot(x, p, "b")
        plt.plot(x, analytical(x, t_end, u, sigma, mu), "k--")

    plt.xlabel(r"Distance [m]")
    plt.ylabel(r"Probability p(x)")
    plt.title(r"Distrubution using FTCS scheme.")
    plt.legend([r"Numerical solution", r"Analytical solution"])
    plt.show()
