#!/usr/bin/python3
import numpy as np

A = np.Matrix([[1, 2, 5], [-1, 3, -1], [2, 1, -2]])  # 3x3 matrix
B = np.Matrix([[1, 2, 5, 2], [-1, 3, -1, -1],
               [2, 1, -2, 1], [1, -1, 1, -1]])  # 4x4 matrix
A = np.linalg.inv(A)
B = np.linalg.inv(B)
print(A)
print(B)
