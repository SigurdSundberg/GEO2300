import numpy as np
from numpy import linalg as LA


A = np.matrix([[1, 0, 5], [-1, 1, 0], [2, 1, -2]])
eigenvalues,  eigenvectors = LA.eig(A)  # v[:,i] = eigenvalue to w[i] A
print(f"Eigenvalues for Matrix A = \n {A}")
print(f"Is as follows: {eigenvalues}")
print(f"and their corresponding eigenvectors are \n {eigenvectors} \n")

"""
Eigenvalues for Matrix A = 
 [[ 1  0  5]
 [-1  1  0]
 [ 2  1 -2]]
Is as follows: [-4.13640529  2.47760887  1.65879642]
and their corresponding eigenvectors are 
 [[-0.69118387 -0.804427   -0.54870228]
 [-0.13456568  0.54441132  0.83288595]
 [ 0.7100401  -0.2377257  -0.07229662]] 
"""
