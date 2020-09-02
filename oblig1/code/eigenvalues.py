import numpy as np
from numpy import linalg as LA


A = np.matrix([[1, 0, 5], [-1, 1, 0], [2, 1, -2]])
eigenvalues,  eigenvectors = LA.eig(A)  # v[:,i] = eigenvalue to w[i] A
print(f"Eigenvalues for Matrix A = \n {A}")
print(f"Is as follows: {eigenvalues}")
print(f"and their corresponding eigenvectors are \n {eigenvectors} \n")
