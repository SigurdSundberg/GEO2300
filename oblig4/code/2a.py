"""
Written by Sigurd Sandvoll Sundberg
Created: 11.11.2020, 15:15 GMT + 1
For the course GEO2300, UiO.
Solutions to problemset 4.
"""
import numpy as np
import matplotlib.pyplot as plt
import warnings

plt.style.use("bmh")
warnings.filterwarnings('error')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "DejaVu Sans",
    "font.sans-serif": ["Helvetica"]})


class Chaos():
    def __init__(self, startingPoint, period, stepLength):
        """Constructor for the chaos system
        Setups our chaos systems initial variables
        Args:
            startingPoint (float): starting point for the 1d case. Can be extended to 3
            period (int): period for the simulation
            stepLength (float): step length for the simulation
        """
        # This program is currently only for a fixed time step
        stepLength = 1
        # This fixes that
        self._n = int(period / stepLength)
        self._period = period
        self._stepLength = stepLength
        self._x = np.zeros(self._n)
        self._x[0] = startingPoint
        self._u0 = startingPoint  # Potential artifact, not sure if needed
        self._timeAtWarning = self._n

    def __call__(self, r):
        """Gives the functionality of the function such that one can F(r) which returns
            wether an overflow is encountered
        Args:
            r (float): varying cases of stability
        Return:
            [boolean]: Wether an overflow is encountered
        """
        if self._x[0] != self._u0:
            self._x = np.zeros(self._n)
            self._x[0] = self._u0
        self._r = r
        for i in range(self._n - 1):
            check = self.advanceOneStep(i)
            if check:
                return True
        return False

    def advanceOneStep(self, index):
        """Advances the time step once

        Args:
            index (int): index of the array
        """
        ix = self._x[index]
        try:  # This treats cases where our equation would be numerically unstable.
            rhs = self._r * ix * (1 - ix)
        except RuntimeWarning:
            print(rf"Overflow in advanceOneStep() for r = {self._r}. Stopping calculations.")
            self._timeAtWarning = index
            return True
        self._x[index + 1] = rhs
        return False

    def _root(self, coefficient=None):
        if coefficient is None:
            coefficient = self._r
        if coefficient == 0:
            return (coefficient, coefficient)
        else:
            root1 = (1 - coefficient) / (2 * coefficient) + ((np.sqrt((1 - coefficient) * (1 - coefficient))) / (2 * coefficient))
            root2 = (1 - coefficient) / (2 * coefficient) - ((np.sqrt((1 - coefficient) * (1 - coefficient))) / (2 * coefficient))
            return (root1, root2)

    def plot(self):
        roots = self._root()
        time = np.linspace(0, self._period, self._n)
        if self._timeAtWarning != self._n:
            plt.plot(time[:self._timeAtWarning], self._x[:self._timeAtWarning], color='k', linestyle='-', markersize=1, marker='o', label=rf"Overflow at iteration: {self._timeAtWarning:d}, for r = {self._r:.1f}", linewidth=.7)

        else:
            plt.plot(time, self._x, color='k', linestyle='-', markersize=1, marker='o', label=rf"r = {self._r:.1f}", linewidth=.7)

        # Plot the roots
        if (roots[0] == roots[1]):
            plt.plot([0, self._timeAtWarning], [roots[0], roots[0]], '--', label=rf"$\lambda_u$ = {roots[0]:.4f}")
        else:
            plt.plot([0, self._timeAtWarning], [roots[0], roots[0]], '--', label=rf"$\lambda^1_u$ = {roots[0]:.4f}")
            plt.plot([0, self._timeAtWarning], [roots[1], roots[1]], '--', label=rf"$\lambda^2_u$ = {roots[1]:.4f}")

    def getArray(self):
        return self._x

    def getRoots(self, r=None):
        return self._root(r)


if __name__ == "__main__":
    # ****************************************
    # Case 0: Initialization for all the other Cases.
    finalTime = 1000
    uInitial = 0.2
    dt = 1
    system = Chaos(uInitial, finalTime, dt)
    rValues = np.arange(0, 10.1, 0.1)
    # ****************************************

    # ****************************************
    # # Case 1: What range of if r \in [0,10.1) do we have a steady state solution? Which root is the steady solution?
    # for i, r in enumerate(rValues):
    #     system(r)
    #     array = system.getArray()
    #     diff = abs(array[-1] - array[-2])
    #     print(diff)
    #     if (diff > 1e-3):
    #         print(rf"Difference between the end points: {diff:.6f}")
    #         break
    # system(rValues[i])
    # system.plot()
    # plt.title("Last value of r which we have a steady state solution")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.legend()
    # # plt.savefig("d.png")
    # plt.show()
    # print(rf"Range is r = [0, ..., {r:.1f}), with precision 1e-16(Machine zero).")
    # input()
    # ****************************************
    print("")
    # ****************************************
    # Case 2: Find an r for which we have a single oscillation?
    i = 30
    for _ in range(len(rValues)):
        r = rValues[i]
        system(r)
        system.plot()
        plt.title("Last value of r which we have a steady state solution")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.show()

        i += 1

    # system(2.7)
    # system.plot()
    # plt.title("A invislbe single oscillation")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.legend()
    # plt.savefig("e1.png")
    # plt.show()
    # system(3.4)
    # system.plot()
    # plt.title("For which we have a single oscillation")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.legend()
    # plt.savefig("e2.png")
    # plt.show()
    # # ****************************************
    # print("")
    # ****************************************
    # # Case 3: Find an r for which we have a double oscillation?
    # system(3.5)
    # system.plot()
    # plt.title("For which we have a double oscillation")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.legend()
    # plt.savefig("f.png")
    # plt.show()
    # # ****************************************
    # print("")
    # # ****************************************
    # # Case 4: For what range of r do we have chaotic oscillations? What are the min and max values of u?
    # system(3.8)
    # system.plot()
    # plt.title("For which we have a chaos")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.legend()
    # plt.savefig("g.png")
    # plt.show()
    # # ****************************************
    # print("")
    # # ****************************************
    # # Case 5: For what r do we have a numerically unstable mapping?
    # system(4.1)
    # system.plot()
    # plt.title("For which we have a chaos")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.legend()
    # plt.savefig("h.png")
    # plt.show()
    # >>> Output: r = 4.1000000000000005
    # ****************************************
