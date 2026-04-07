#!/usr/bin/env python3
"""Verify the C steepest-descent output against NumPy/SciPy."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

try:
    from scipy.linalg import solve as scipy_solve
except Exception:
    scipy_solve = None


def load_vector(path: Path) -> np.ndarray:
    return np.loadtxt(path, dtype=np.float64)


def load_matrix(path: Path) -> np.ndarray:
    return np.loadtxt(path, dtype=np.float64)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", nargs="?", default="verify")
    args = parser.parse_args()

    prefix = Path(args.prefix)
    a_path = prefix.with_name(prefix.name + "_A.txt")
    b_path = prefix.with_name(prefix.name + "_b.txt")
    x_path = prefix.with_name(prefix.name + "_x.txt")
    x_true_path = prefix.with_name(prefix.name + "_x_true.txt")

    A = load_matrix(a_path)
    b = load_vector(b_path)
    x_c = load_vector(x_path)

    x_ref = None
    if x_true_path.exists():
        x_ref = load_vector(x_true_path)
    elif scipy_solve is not None:
        x_ref = scipy_solve(A, b)
    else:
        x_ref = np.linalg.solve(A, b)

    residual = np.linalg.norm(A @ x_c - b)
    rel_err_ref = np.linalg.norm(x_c - x_ref) / np.linalg.norm(x_ref)

    print(f"matrix shape      : {A.shape}")
    print(f"residual norm     : {residual:.6e}")
    print(f"relative error    : {rel_err_ref:.6e}")
    print(f"reference solver  : {'x_true file' if x_true_path.exists() else 'SciPy/NumPy solve'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
