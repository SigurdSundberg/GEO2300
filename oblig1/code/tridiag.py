import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA


def tridiag_solver_Gauss(n):
    """solve the problem using gaussian elim

    Args:
        n (int): number of points

    Returns:
        1d-array: solution to the problem
    """
    d = np.zeros(n)  # diagonal
    d.fill(-2)  # fill the array
    solution = np.zeros(n)  # vector x
    d[0] = d[n-1] = 1  # first element and last element

    b = vec_b(n)

    # Forwards sub
    for i in range(2, n-1):
        d[i] = -2 - 1/d[i-1]
        b[i] = b[i] - b[i-1]/d[i-1]
    # backwards sub
    solution[n-2] = b[n-2]/d[n-2]
    for i in range(n-3, 0, -1):
        solution[i] = (b[i] - solution[i+1])/d[i]
    return solution


def matrix_solution(n):
    """solve the problem with matrix multiplication and inversion

    Args:
        n (int): number of points

    Returns:
        1d-array: solution to the answer
    """
    A = np.zeros([n, n])
    A[0, 0] = A[n-1, n-1] = 1
    solution = -1
    for i in range(1, n-1):
        A[i, i+1] = 1
        A[i, i] = -2
        A[i, i-1] = 1
    A_inv = LA.inv(A)
    solution = A_inv@vec_b(n)
    return solution


def vec_b(n):
    """Generate vector b. b[0] = b[n-1] = 0

    Args:
        n (int): number of points

    Returns:
        1d-array: returns the array with values in b
    """
    b = np.zeros(n)  # empty vector
    delta_y = h/(n-1)  # steplength
    b.fill(dpdx / my * delta_y * delta_y)
    b[0] = b[n-1] = 0
    # print(f"b = {b}")
    return b


def vec_x(n):
    """set up x vector

    Args:
        n (int): number of points

    Returns:
        1d-array: vector of x-values
    """
    x = np.linspace(0, h, n)
    # print(f"x = {x}")
    return x


def plotter(x, data, solver):
    """Plotting the data

    Args:
        x (1d-array): x-values
        data (1d-array): data values
        solver (functions): which solver for the data set is used
    """
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica"]})

    plt.figure()
    plt.plot(x, data)
    # plt.plot(x, analytical_solution(x))
    plt.xlabel(r"Velocity $\left[\displaystyle\frac{m}{s}\right]$")
    plt.ylabel(r"Height $[\displaystyle m]$")
    plt.title(solver)


def analytical_solution(x):
    """finds the analytical solution to the problem

    Args:
        x (1d-array): x-values

    Returns:
        1d-array: solution to the problem
    """
    return -1


if __name__ == "__main__":
    list_of_n = (5, 100)
    solver = ("Gauss elim", "matrix")
    dpdx = -200  # Pa/m
    h = 10  # m
    my = 300  # Pa*s
    for ele in list_of_n:
        print(ele)
        x = vec_x(ele)
        data1 = tridiag_solver_Gauss(ele)
        data2 = matrix_solution(ele)
        plotter(data1, x, solver[0])
        plotter(data2, x, solver[1])
    plt.show()
