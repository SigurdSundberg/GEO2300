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


def Diffusion(x, t, u, sigma, mu, c, D):
    if t == 0:
        leading_term = 1/(sigma * np.sqrt(2*np.pi))
        divisor = (2*sigma*sigma)
        return leading_term * np.exp(-(x-u*t-mu)*(x-u*t-mu) / divisor)
    else:
        return (1/(2*np.sqrt(np.pi*D*t)))*np.exp(-(x-u*t-mu)*(x-u*t-mu)/(4*D*t))


def No_Diffusion(x, t, u, sigma, mu, c, D):
    leading_term = 1/(sigma * np.sqrt(2*np.pi))
    divisor = (2*sigma*sigma)
    return leading_term * np.exp(-(x-u*t-mu)*(x-u*t-mu) / divisor)


def FTCS_matrix(n, c, s):
    A = np.zeros([n, n])
    diagonal = 1 - 2*s
    upper = s - c/2
    lower = s + c/2
    for i in range(1, n-1):
        A[i, i] = diagonal
        A[i, i-1] = lower
        A[i, i+1] = upper
    A[0, 0] = diagonal
    A[-1, -1] = diagonal
    A[0, 1] = upper
    A[-1, -2] = lower
    return A


X = np.linspace(0, 1000, 1000)

end_times = [0, 500, 1000, 1500]  # [s]
# end_times = [500]  # [s]
C = 0.1  # courant number
mu = 100  # [m] Relative mean
sigma = 10  # [m] Standard deviation
D = 0.5  # [m/s]

dx = 5
n = ceil(1000/dx)
x = np.linspace(0, 1000, n+1)  # meters
u = 0.5  # m/s

dt = C/u * dx
s = D * dt / (dx*dx)

A = FTCS_matrix(n+1, C, s)

pn = u * sigma / D
print(f"Peclet number: {pn:.3f}")

for t_end in end_times:
    t = 0
    p = pdf(x, sigma, mu)
    while t < t_end:
        t += dt
        p = np.dot(A, p)
    if t_end == 0:
        plt.plot(X, No_Diffusion(X, t_end, u, sigma, mu, C, D), "b")
    else:
        plt.plot(X, Diffusion(X, t_end, u, sigma, mu, C, D), "k")
        plt.plot(X, No_Diffusion(X, t_end, u, sigma, mu, C, D), "r")
        plt.plot(x, p, "g--", linewidth=2)
plt.xlabel(r"Distance [m]")
plt.ylabel(r"Probability p(x)")
plt.title(r"Distrubution including diffusion.")
plt.legend([r"Initial distrubution.", r"Distrubution with diffusion.",
            r"Distrubution with no diffusion.", r"Numerical distrubution including diffusion."])
plt.show()
