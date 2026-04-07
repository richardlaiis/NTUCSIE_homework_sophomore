# solution_checker.py
import numpy as np
from scipy.linalg import solve

A = np.loadtxt('A', dtype=np.float64)
b = np.loadtxt('b', dtype=np.float64)

x = solve(A, b)
print(x)