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
    d = np.zeros(n)  # diagonal. Could be done with np.full also.
    d.fill(-2)  # fill the array
    solution = np.zeros(n)  # vector v
    d[0] = d[n-1] = 1  # Set first and last element of diag

    b = vec_b(n)  # Setup vector b

    # Forwards sub
    for i in range(2, n-1):
        d[i] = -2 - 1/d[i-1]
        b[i] = b[i] - b[i-1]/d[i-1]
    # backwards sub
    solution[n-2] = b[n-2]/d[n-2]  # End points
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
    A = np.zeros([n, n])  # Matrix setup
    A[0, 0] = A[n-1, n-1] = 1
    for i in range(1, n-1):
        A[i, i+1] = 1
        A[i, i] = -2
        A[i, i-1] = 1
    A_inv = LA.inv(A)  # Inversion
    return A_inv@vec_b(n)  # Matrix multiplication


def vec_b(n):
    """Generate vector b. b[0] = b[n-1] = 0

    Args:
        n (int): number of points

    Returns:
        1d-array: returns the array with values in b
    """
    delta_y = h/(n-1)  # steplength
    b = np.full(n, dpdx / my * delta_y * delta_y)  # Setup
    b[0] = b[n-1] = 0  # endpoints
    return b


def vec_x(n):
    """set up x vector

    Args:
        n (int): number of points

    Returns:
        1d-array: vector of x-values
    """
    x = np.linspace(0, h, n)
    return x


def plotter(x, data, solver):
    """Plotting the data

    Args:
        x (1d-array): x-values
        data (1d-array): data values
        solver (functions): which solver for the data set is used
    """
    plt.figure()
    plt.plot(data, x, label=(solver + " n = " + str(ele)))
    plt.title("Solution from " + solver + " method")
    plot_analytical_solution()
    plt.xlabel(r"Velocity $\left[\displaystyle\frac{m}{s}\right]$")
    plt.ylabel(r"Height $[\displaystyle m]$")
    plt.legend()
    plt.savefig("./plot/plot_" + solver + str(ele))
    plt.title(solver)


def plot_analytical_solution():
    """Plots the analytical solution and returns the values

    Returns:
        1d-array: return values of the solution on vector form. 
    """
    const = dpdx / (2*my)
    n = 1000
    x = np.linspace(0, h, n)
    solution = const * (x*x-h*x)
    plt.plot(solution, x, "k--", label="Analytical")
    return solution


if __name__ == "__main__":
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "DejaVu Sans",
        "font.sans-serif": ["Helvetica"]})

    # To solve for both N, and possibly more values if needed.
    list_of_n = (5, 100)
    # to get two plots, using either method for solving
    solver = ("Gauss Elimination", "Matrix equation")
    dpdx = -200  # Pa/m
    h = 10  # m
    my = 300  # Pa*s
    for ele in list_of_n:
        x = vec_x(ele)
        data1 = tridiag_solver_Gauss(ele)
        data2 = matrix_solution(ele)
        plotter(x, data1, solver[0])
        plotter(x, data2, solver[1])
        print(
            f"The maximum speed using numerical solution for n = {ele} is {max(data1):5.3f}[m/s]")
        print(
            f"The maximum speed using the analytical solution for n = 1000 is {max(plot_analytical_solution()):5.3f}[m/s]")
    plt.show()

"""
The maximum speed using numerical solution for n = 5 is 8.333[m/s]
The maximum speed using the analytical solution for n = 1000 is 8.333[m/s]
The maximum speed using numerical solution for n = 100 is 8.332[m/s]
The maximum speed using the analytical solution for n = 1000 is 8.333[m/s]

"""
