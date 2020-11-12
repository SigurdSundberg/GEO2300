import numpy as np
import matplotlib.pyplot as plt
from scipy import special

if __name__ == "__main__":
    cms_to_ms = 1/(100*100)  # cm^2/s to m^2/s
    day_to_s = 24 * 60 * 60  # days to seconds

    def diffusion(t, x, c):
        D = 1 * cms_to_ms  # cm^2/s
        return c*special.erfc(x/(2*np.sqrt(D*t)))

    c0 = 0.5
    t_list = np.array([1, 10, 100]) * day_to_s
    x = np.linspace(0, 100, 2000)

    x_list = np.array([100, 500, 1000, 1500, 2000])
    dt = 1
    c = 0.5
    c1 = 0
    t_old = []
    for x_value in x_list:
        c1 = 0
        x_array = np.linspace(0, x_value, 5000)
        t_start = 0
        while c1 < 0.25:
            t_start += dt
            c1 = diffusion(t_start*day_to_s, x_array, c0)
            c1 = c1[-1]
        print(
            f"Number of days needed for the concentration to reach 25% is: {t_start:d} days, with x = {x_value:d}")
        t_old.append(t_start)
    for t in t_list:
        plt.plot(x, diffusion(t, x, c0),
                 label=f"Diffusion for t = {t/day_to_s:.0f}")
    plt.savefig("FTCS2.pdf")
    plt.rcParams.update({
        "text.usetex": True,
        "font.sans-serif": ["Helvetica"],
        "font.family": "DejaVu Sans"})
    plt.title(r"Analytical solution for Diffusion equation")
    plt.xlabel(r"Distance x $\left[ m\right]$")
    plt.ylabel(r"Concentration $[-]$")
    plt.legend()
    plt.savefig("FTCS2.pdf")
    plt.show()

    q = np.log10(x_list)
    p = np.log10(t_old)
    m, b = np.polyfit(q, p, 1)
    plt.plot(q, p, "r")
    plt.plot(q, p, "ko", label=r"datapoints")
    plt.annotate(fr"The incline of the line goes as m = {m:.2f}", [2, 5])
    plt.title(r"Log/Log plot of distance vs time")
    plt.xlabel(r"$\log_{10} (x)$")
    plt.ylabel(r"$\log_{10} (t)$")
    plt.legend()
    plt.savefig("Relation.pdf")
    plt.show()
"""
Number of days needed for the concentration to reach 25% is: 1273 days
Number of days needed for the concentration to reach 25% is: 31802 days
Number of days needed for the concentration to reach 25% is: 127206 days
Number of days needed for the concentration to reach 25% is: 286213 days
Number of days needed for the concentration to reach 25% is: 508822 days
"""
