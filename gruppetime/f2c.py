import numpy as np
import matplotlib.pyplot as plt
import datetime as dt


def f2c(f):
    return 5 / 9 * (f - 32)


def test_f2c():
    f = 32
    c = 0
    out = f2c(f) == c
    msg = "Output is {} should be {}".format(out, c)
    assert out, msg


test_f2c()

data = np.loadtxt("dallas_temp.txt")
month = data[:, 0]
day = data[:, 1]
year = data[:, 2]
temp_f = data[:, 3]

temp_c = f2c(temp_f)

dates = [dt.date(int(year[i]), int(month[i]), int(day[i])) for i in range(len(temp_f))]
fig, ax = plt.subplots()
ax.plot_date(dates, temp_c, marker="+", markersize=4)
# time = np.linspace(0, len(temp), len(temp))
plt.xlabel("time since first measurment [d]")
plt.ylabel("temperature in c [C]")
plt.show()

for i, k in enumerate(temp_f):
    if k == -99:
        print(i)
        temp_f[i] = (temp_f[i - 1] + temp_f[i + 1]) / 2
fig, ax2 = plt.subplots()
ax2.plot_date(dates, f2c(temp_f), marker="+", markersize=4)
plt.xlabel("time since first measurment [d]")
plt.ylabel("temperature in c [C]")
plt.show()
