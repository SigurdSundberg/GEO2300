import numpy as np
import matplotlib.pyplot as plt
import time

print("This could be very well be an implicit scheme that has been implemented instead of an explicit scheme as it should be")
print("This is most an implicit scheme. HOW THE FUCK DID I MANAGE THAT?")


def B(x):
    return 1.0 * np.exp(-(x*x)/(50*50))


def analytical_solution(t):
    A = 1.0
    std = 50
    c = 700
    x = np.linspace(-1400, 300, 400)
    X = x + c * t
    y = A * np.exp(-(X*X)/(std*std))
    return x, y


def tridiag_solver_Gauss():
    c_0 = 700  # km/h

    N = 400
    C = 0.01
    n = N - 2

    d = 1
    L = -C/2
    U = C/2

    d0 = (1-C/2)
    dn = (1+C/2)

    X = np.linspace(-1400, 300, N)
    dx = X[1]-X[0]
    X = X[1:-1]

    dt = C*dx/c_0

    # [1:-1]
    b = B(X)
    print(b)
    tL = [0, 0.5, 1, 1.5, 2]
    for t_end in tL:
        t = 0
        x, y = analytical_solution(t_end)
        b = B(X)
        print(b)
        b_sol = np.zeros(N-2)
        stime = time.time()
        while t < t_end:
            b_sol[0] = U * b[1] + b[0] * d0
            for i in range(1, n-1):
                b_sol[i] = d * b[i] + U * (b[i+1] - b[i-1])
            b_sol[n-1] = L * b[n-2] + b[n-1] * dn
            b = b_sol
            t += dt
        etime = time.time()
        print(b.argmax())
        # print(etime - stime)
        plt.plot(X, b, 'r-.')
        plt.plot(x, y, 'k--')
    plt.savefig("stable0.01.pdf")
    plt.show()

    # return solution
tridiag_solver_Gauss()
