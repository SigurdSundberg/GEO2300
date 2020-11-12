import numpy as np
import matplotlib.pyplot as plt


def x(p):
    return 2*p + 1


n = 10*100*1000*10

X = []
Y = []
for i in range(1, n):
    Y.append(x(i))
    X.append(i)

x = np.array(X)
y = np.array(Y)

log_x = np.log(x)
log_y = np.log(y)
np.insert(log_x, 0, 0)
np.insert(log_y, 0, 0)
curve_fit = np.polyfit(log_x, log_y, 1)
print(curve_fit)
"""  
[0.99999411 0.69323699]
>>> 0.99999411*exp(0.69323699*x)
"""
