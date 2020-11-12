import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA


def H(x):
    return 1.0 * np.exp(-(x*x)/(50*50))


def initialize_tridiag(diagonal, offDiagonalL, offDiagonalU, n):
    A = np.zeros([n, n])

    for i in range(1, n-1):
        A[i, i] = diagonal
        A[i, i-1] = offDiagonalL
        A[i, i+1] = offDiagonalU
    return A


def analytical_solution(t):
    A = 1.0
    std = 50
    c = 700
    x = np.linspace(-1400, 300, 400)
    X = x + c * t
    y = A * np.exp(-(X*X)/(std*std))
    return x, y


def solver(t_end):
    c_0 = 700

    N = 400
    C = 0.01
    n = N - 2

    d = 1
    nd = C/4
    A = initialize_tridiag(d, -nd, nd, n)
    B = initialize_tridiag(d, nd, -nd, n)
    X = np.linspace(-1400, 300, N)[1:-1]

    dx = X[1]-X[0]
    dt = C*dx/c_0

    # Apply boundaries
    A[0, 0] = (d-nd)
    A[n-1, n-1] = d+nd
    A[0, 1] = nd
    A[n-1, n-2] = -nd
    B[0, 0] = d + nd
    B[n-1, n-1] = d - nd
    B[0, 1] = - nd
    B[n-1, n-2] = nd
    Binv = LA.inv(B)

    t = 0
    x, y = analytical_solution(t_end)
    h = H(X)
    while t <= t_end:
        t += dt
        h = np.dot(Binv, np.dot(A, h))
    return [(X, h), (x, y)]


if __name__ == "__main__":
    plt.rcParams.update({
        "text.usetex": True,
        "font.sans-serif": ["Helvetica"],
        "font.family": "DejaVu Sans"})
    plt.title(r"Wave equation using CN")
    plt.xlabel(r"Position $\left[ m\right]$")
    plt.ylabel(r"Amplitude $[-]$")
    color = ["c", "r", "g", "b", "m"]
    i = 0
    for t_end in [0, 0.5, 1, 1.5, 2]:
        [X, h], [x, y] = solver(t_end)
        plt.plot(
            X, h, color[i], label=fr"Numerical solution for t = {t_end:.1f}", linewidth=1)
        plt.plot(
            x, y, color[i] + '--', label=fr"Analytical solution for t = {t_end:.1f}", linewidth=1)
        i += 1

    # plt.legend()
    plt.savefig("CN.pdf")
    plt.show()
