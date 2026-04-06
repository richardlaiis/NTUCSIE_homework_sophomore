import numpy as np
import scipy.linalg
import sys

def verify_steepest_descent(d):
    # Generate the same SPD matrix as the C code
    # A_ij = 1.0 if i != j else d + 1.0
    A = np.ones((d, d)) + np.eye(d) * d
    
    # b_i = i + 1.0
    b = np.arange(1, d + 1, dtype=float)
    
    # Solve Ax = b
    x = scipy.linalg.solve(A, b)
    
    print(f"Dimension: {d} x {d}")
    print("SciPy Solution for first 5 elements:")
    for i in range(min(5, d)):
        print(f"x[{i}] = {x[i]:.6f}")

if __name__ == "__main__":
    d = 100
    if len(sys.argv) > 1:
        d = int(sys.argv[1])
    verify_steepest_descent(d)
