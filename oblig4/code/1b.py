"""
Written by Sigurd Sandvoll Sundberg
Created: 09.11.2020, 16:00 GMT + 1
For the course GEO2300, UiO.
Solutions to problemset 4.

Solves Burgers equation using the FTCS scheme.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import warnings
plt.style.use("bmh")
warnings.filterwarnings('error')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "DejaVu Sans",
    "font.sans-serif": ["Helvetica"]})


class Analytical():
    def __init__(self, a, v):
        # Constructor
        self._v = v
        self._a = a

    def __call__(self, x, t):
        # Function call
        a = self._a
        v = self._v
        return .5 * a - .5 * a * np.tanh(a * ((x - (.5 * a * t)) / (4 * v)))

    def plot(self):
        pass


class Burgers():
    def __init__(self, s, c, n):
        self._s = s  # Problem defined variable
        self._c = c  # Problem defined variable
        self._n = n  # Dimensionality
        self._A = np.zeros([n, n])  # Matrix A
        self._B = np.zeros([n, n])  # Matrix B

    def dirichlet(self, Boundaries):
        self.bound = np.zeros(self._n)
        self.bound[0] = Boundaries[0] * s + Boundaries[0] * c
        self.bound[-1] = Boundaries[1] * s + Boundaries[1] * c

    def neumann(self):
        raise NotImplementedError

    def setup(self, Boundaries, Type):
        if Type == "D":
            self.dirichlet(Boundaries)
        elif Type == "N":
            self.neumann()
        else:
            sys.exit("No legal boundaries given")
        n = self._n
        c = self._c
        s = self._s
        for i in range(1, n - 1):
            self._A[i, i - 1] = c
            self._A[i, i + 1] = -c
            self._B[i, i - 1] = s
            self._B[i, i] = 1 - 2 * s
            self._B[i, i + 1] = s
        # Matrix boundaries
        self._A[0, 1] = -c
        self._B[0, 0] = 1 - 2 * s
        self._B[0, 1] = s
        self._A[-1, -2] = c
        self._B[-1, -1] = 1 - 2 * s
        self._B[-1, -2] = s

    def FTCS(self, dt, t_initial, t_end, U):
        self._T = t_end
        self._dt = dt
        self._u = U  # Initial function

        def Evec(): return 0.5 * (self._u * self._u)

        t = t_initial
        E = Evec()
        while t < self._T:
            t += dt
            self._u = np.dot(self._A, E) + \
                np.dot(self._B, self._u) + self.bound
            try:
                E = Evec()
            except RuntimeWarning:
                print("Overflow in E, returning to __main__")
                break

    def plot(self, x, nu=None, a=None):
        plt.plot(x, self._u, '--', label=rf"Numerical at t = {self._T:d}")
        plt.title(
            rf"$a = {a:.1f}, \nu = {nu:.2f}, s = {self._s:.2f}$")
        plt.xlabel(r"Distance [km]")
        plt.ylabel(r"Height [-]")


if __name__ == '__main__':

    # ****************************************************** #
    # Case 1: a = 2, v = 1, s = 0.4
    a, v = (2, 1)  # Amplitude and viscocity
    s = 0.4        # Variable defined by properties of the system
    B = (2, 0)     # Boundaries (L, U)
    print(f"Case 1; a = {a}, v = {v}, s = {s}")

    # Domain setup
    n = 750
    time_initial = 1
    t_list = [10, 20]
    x = np.linspace(-10, 30, n)
    dx = x[1] - x[0]

    # Fixed variables
    dt = (s * dx * dx) / v
    c = dt / (2 * dx)

    start_time = time.time()
    wave = Analytical(a, v)
    initial_wave = wave(x, time_initial)

    solver = Burgers(s, c, n)
    solver.setup(B, "D")
    print("Setup --- %s seconds ---" % (time.time() - start_time))
    plt.plot(x, initial_wave, 'b', label=rf"Initial wave")
    for t in t_list:
        start_time = time.time()
        plt.plot(x, wave(x, t), 'k', label=rf"Analytical at t = {t:d}", linewidth=2)
        solver.FTCS(dt, time_initial, t, initial_wave)
        solver.plot(x, v, a)
        print("t = %d --- %s seconds ---" %
              (t, (time.time() - start_time)))
    plt.legend()
    plt.savefig("fig1.png")
    plt.show()
    # ****************************************************** #
    print("")
    # ****************************************************** #
    # Case 2: a = 2, v = 0.1, s = 0.4
    a, v = (2, 0.1)  # Amplitude and viscocity
    s = 0.4          # Variable defined by properties of the system
    B = (2, 0)       # Boundaries (L, U)
    print(f"Case 2; a = {a}, v = {v}, s = {s}")

    # Domain setup
    n = 1500
    t_list = [10, 20]
    x = np.linspace(-10, 30, n)
    dx = x[1] - x[0]

    # Fixed variables
    dt = (s * dx * dx) / v
    c = dt / (2 * dx)

    start_time = time.time()
    wave = Analytical(a, v)
    initial_wave = wave(x, time_initial)

    solver = Burgers(s, c, n)
    solver.setup(B, "D")
    print("Setup --- %s seconds ---" % (time.time() - start_time))
    plt.plot(x, initial_wave, label=rf"Initial wave")
    for t in t_list:
        start_time = time.time()
        plt.plot(x, wave(x, t), 'k', label=rf"Analytical at t = {t:d}", linewidth=2)
        solver.FTCS(dt, time_initial, t, initial_wave)
        solver.plot(x, v, a)
        print("t = %d --- %s seconds ---" %
              (t, (time.time() - start_time)))
    plt.legend()
    plt.savefig("fig2.png")
    plt.show()
    # ****************************************************** #
    print("")
    # ****************************************************** #
    # Case 3: a = 2, v = 0.1, s = 0.6
    a, v = (2, 0.1)  # Amplitude and viscocity
    s = 0.6         # Variable defined by properties of the system
    B = (2, 0)       # Boundaries (L, U)
    print(f"Case 3; a = {a}, v = {v}, s = {s}")

    # Domain setup
    n = 1500
    t_list = [10, 20]
    x = np.linspace(-10, 30, n)
    dx = x[1] - x[0]

    # Fixed variables
    dt = (s * dx * dx) / v
    c = dt / (2 * dx)

    start_time = time.time()
    wave = Analytical(a, v)
    initial_wave = wave(x, time_initial)

    solver = Burgers(s, c, n)
    solver.setup(B, "D")
    print("Setup --- %s seconds ---" % (time.time() - start_time))
    plt.plot(x, initial_wave, label=rf"Initial wave")
    for t in t_list:
        start_time = time.time()
        plt.plot(x, wave(x, t), 'k', label=rf"Analytical at t = {t:d}", linewidth=2)
        solver.FTCS(dt, time_initial, t, initial_wave)
        solver.plot(x, v, a)
        print("t = %d --- %s seconds ---" %
              (t, (time.time() - start_time)))
    plt.legend()
    plt.savefig("fig3.png")
    plt.show()
    # ****************************************************** #
