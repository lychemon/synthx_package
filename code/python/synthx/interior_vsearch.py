"""
Compute cell weights and test diagonal attainability for interior points.
"""

from __future__ import annotations

import itertools
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import lu_factor, lu_solve

LP = dict(method="highs")
TOL = 1e-9


def cell_factor(A, X0):
    """
    Factor the cell matrix once for repeated barycentric solves.
    """
    M = np.vstack([X0[:, A], np.ones(len(A))])
    return lu_factor(M)


def barycentric(A, X0, X1, luf=None):
    if luf is None:
        M = np.vstack([X0[:, A], np.ones(len(A))])
        return np.linalg.solve(M, np.append(X1, 1.0))
    return lu_solve(luf, np.append(X1, 1.0))


def admissible(A, X0, X1, luf=None):
    w = barycentric(A, X0, X1, luf)
    return (w > TOL).all(), w


def rho_rows(A, X0, luf=None):
    """
    Return the diagonal cell inequalities for donors outside A.
    """
    if luf is None:
        luf = cell_factor(A, X0)
    K, J = X0.shape
    rows = []
    for j in range(J):
        if j in A:
            continue
        beta = lu_solve(luf, np.append(X0[:, j], 1.0))
        diagQ = X0[:, j] ** 2 - sum(beta[i] * X0[:, A[i]] ** 2 for i in range(len(A)))
        rows.append(diagQ)
    return np.array(rows) if rows else np.zeros((0, K))


def theta_diag(A, X0, luf=None):
    """
    Maximize the common slack in the diagonal cell inequalities.
    """
    R = rho_rows(A, X0, luf)
    K = X0.shape[0]
    if R.shape[0] == 0:
        return np.inf
    c = np.zeros(K + 1)
    c[-1] = -1.0                                   # maximize theta
    A_ub = np.hstack([-R, np.ones((R.shape[0], 1))])   # -rho_j'd + theta <= 0
    b_ub = np.zeros(R.shape[0])
    A_eq = np.zeros((1, K + 1))
    A_eq[0, :K] = 1.0
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0],
                  bounds=[(0, None)] * K + [(None, None)], **LP)
    return -res.fun if res.status == 0 else -np.inf
