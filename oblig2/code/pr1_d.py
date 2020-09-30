import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import time


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

    diagonal = 1
    offDiagonalL = -C/2
    offDiagonalU = C/2
    A = initialize_tridiag(diagonal, offDiagonalL, offDiagonalU, n)
    X = np.linspace(-1400, 300, N)[1:-1]

    dx = X[1]-X[0]
    dt = C*dx/c_0
    # Apply boundaries
    A[0, 0] = (1-C/2)
    A[n-1, n-1] = 1+C/2

    A[0, 1] = C/2
    A[n-1, n-2] = - C/2
    t = 0
    x, y = analytical_solution(t_end)
    h = H(X)
    print(h)
    stime = time.time()
    while t <= t_end:
        t += dt
        h = np.dot(A, h)
    etime = time.time()
    # print(etime - stime)
    return [(X, h), (x, y)]


for t_end in [0, 0.5, 1, 1.5, 2]:
    print(t_end)
    [X, h], [x, y] = solver(t_end)
    plt.plot(X, h, 'r-.')
    plt.plot(x, y, 'k--')
# plt.savefig("FTCS.pdf")
plt.show()
