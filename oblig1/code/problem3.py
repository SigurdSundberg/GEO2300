import sympy as sp

s = sp.Symbol('s')
A = sp.Matrix([[1+2*s, -s, 0], [-s, 1+2*s, -s], [0, -s, 1+2*s]])
A_inv = A.inv()
print("A inverse inserted for s = 0.5 gives")
sp.pprint(A_inv.evalf(subs={s: 0.5}))
"""
A inverse inserted for s = 0.5 gives
|0.535714285714286   0.142857142857143  0.0357142857142857|
|                                                         |
|0.142857142857143   0.571428571428571  0.142857142857143 |
|                                                         |
|0.0357142857142857  0.142857142857143  0.535714285714286 |
"""
