"""
Written by Sigurd Sandvoll Sundberg
Created: 09.11.2020, 15:30 GMT + 1
For the course GEO2300, UiO. 
Solutions to problemset 4. 
"""
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("bmh")
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
        return a / 2 - a / 2 * np.tanh(a * ((x - (a * t / 2)) / (4 * v)))


def plotting():
    # Fixed variables
    t_list = [1, 10, 20]
    x = np.linspace(-10, 30, 500)

    # ****************************************************** #
    # Case 1: a = 2, v = 1.
    a, v = (2, 1)
    f = Analytical(a, v)
    for t in t_list:
        plt.plot(x, f(x, t), label=rf"t = {t:d}")
    plt.title(r"$a = 2, \nu = 1$")
    plt.xlabel(r"Distance [km]")
    plt.ylabel(r"Height [m]")
    plt.legend()
    plt.savefig("1a1.png")
    plt.show()
    # ****************************************************** #

    # ****************************************************** #
    # Case 2: a = 2, v = 0.1.
    a, v = (2, 0.1)
    f = Analytical(a, v)
    for t in t_list:
        plt.plot(x, f(x, t), label=rf"t = {t:d}")
    plt.title(r"$a = 2, \nu = 0.1$")
    plt.xlabel(r"Distance [km]")
    plt.ylabel(r"Height [m]")
    plt.legend()
    plt.savefig("1a2.png")
    plt.show()
    # ****************************************************** #

    # ****************************************************** #
    # Case 3: a = 2, v = 0.1 and v = 1
    # c = ['b', 'r', 'k']
    # i = 0
    # a, v1 = (2, 0.1)
    # a, v2 = (2, 1)
    # f1 = Analytical(a, v1)
    # f2 = Analytical(a, v2)
    # for t in t_list:
    #     plt.plot(x, f1(x, t), c[i] + '--', label=rf"$\nu = 0.1$")
    #     plt.plot(x, f2(x, t), c[i], label=rf"$\nu = 1$")
    #     i += 1
    # plt.title(r"$a = 2, \nu = 0.1 \land \nu = 1$")
    # plt.xlabel(r"Distance [km]")
    # plt.ylabel(r"Height [m]")
    # plt.legend()
    # plt.savefig("1a3.png")
    # plt.show()
    # ****************************************************** #

    # ****************************************************** #
    # Case 4: a = 4, v = 0.1.
    a, v = (4, 0.1)
    f = Analytical(a, v)
    for t in t_list:
        plt.plot(x, f(x, t), label=rf"t = {t:d}")
    plt.title(r"$a = 4, \nu = 0.1$")
    plt.xlabel(r"Distance [km]")
    plt.ylabel(r"Height [m]")
    plt.legend()
    plt.savefig("1a4.png")
    plt.show()
    # ****************************************************** #


if __name__ == '__main__':
    print("Plotting 1a.")
    plotting()
else:
    print("Plotting from the part 1a) of problemset 4 GEO2300")
    plotting()
