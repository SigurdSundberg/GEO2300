import numpy as np
import matplotlib.pyplot as plt
from scipy import special

if __name__ == "__main__":
    cms_to_ms = (60*60*24)/(100*100)  # cm^2/s to m^2/s
    day_to_s = 1  # 24 * 60 * 60  # days to seconds

    def diffusion(t, x, c):
        D = 1 * cms_to_ms  # cm^2/s
        return c*special.erfc(x/(2*np.sqrt(D*t)))

    N = 100
    n = N
    xmin = 0
    xmax = 100
    X = np.linspace(xmin, xmax, N*5)
    dx = (X[1]-X[0])*5
    x = np.linspace(xmin, xmax, n)
    c0 = 0.5
    D = 1 * cms_to_ms
    s = 0.4
    dt = s/D * dx*dx
    A = np.zeros([n, n])
    mO = s
    mD = (1-2*s)
    C = np.zeros(n)
    cU = 0
    cL = 0.5
    C[0] = cL

    # Setup A
    A[0, 0] = 1
    A[-1, -1] = 1
    for i in range(1, n-1):
        A[i, i] = mD
        A[i, i-1] = mO
        A[i, i+1] = mO

    t_list = np.array([1, 10, 100])
    plt.savefig("FTCS2.pdf")
    plt.rcParams.update({
        "text.usetex": True,
        "font.sans-serif": ["Helvetica"],
        "font.family": "DejaVu Sans"})
    plt.title(r"Analytical solution for Diffusion equation")
    plt.xlabel(r"Distance x $\left[ m\right]$")
    plt.ylabel(r"Concentration $[-]$")
    for t in t_list:
        ts = 0
        while ts <= t:
            ts += dt
            C = np.dot(A, C)
        plt.plot(x, C, "--", label=f"Numerical solution for t = {t:.0f}")
        plt.plot(X, diffusion(t*day_to_s, X, c0), 'k-.',
                 label=f"Analytical solution for t = {t:.0f}")
    plt.legend()
    # plt.savefig("exNUMFTCS.pdf")
    plt.show()
