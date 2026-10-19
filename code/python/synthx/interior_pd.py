"""
Find the best interior cell under a full positive definite metric.
"""

from __future__ import annotations

from collections import deque
import numpy as np
import cvxpy as cp
from scipy.linalg import lu_solve

from .interior_vsearch import admissible, barycentric, theta_diag, cell_factor
from .interior_cellwalk import seed_cell, flip_neighbor, cell_neighbors

TOL = 1e-9


def Q_matrices(A, X0, luf=None):
    """
    Return Q_j = X_j X_j.T - sum_i beta_i X_i X_i.T for outside donors.
    """
    A = list(A)
    if luf is None:
        luf = cell_factor(A, X0)
    K, J = X0.shape
    out = []
    for j in range(J):
        if j in A:
            continue
        beta = lu_solve(luf, np.append(X0[:, j], 1.0))
        Qj = np.outer(X0[:, j], X0[:, j]) - sum(
            beta[i] * np.outer(X0[:, A[i]], X0[:, A[i]]) for i in range(len(A)))
        out.append(0.5 * (Qj + Qj.T))                 # symmetrize against round-off
    return out


def theta_pd(A, X0, luf=None, maxit=200, tol=1e-9):
    """
    Maximize the common slack in the PD cell inequalities.
    """
    Qs = Q_matrices(A, X0, luf)
    if not Qs:
        return np.inf
    K = X0.shape[0]
    V = cp.Variable((K, K), PSD=True)
    slack = cp.Variable()
    constraints = [cp.trace(V) == K]
    constraints += [cp.trace(Q @ V) >= slack for Q in Qs]
    problem = cp.Problem(cp.Maximize(slack), constraints)
    problem.solve(solver=cp.CLARABEL, max_iter=maxit,
                  tol_gap_abs=tol, tol_gap_rel=tol, tol_feas=tol)
    if problem.status != cp.OPTIMAL:
        raise RuntimeError("PD cell SDP failed: " + problem.status)
    return float(problem.value)


def pd_attainable(A, X0):
    """
    Test diagonal attainability first, then allow a full PD metric.
    """
    if theta_diag(list(A), X0) > TOL:
        return True
    return theta_pd(list(A), X0) > TOL


def walk_attainable_cells_pd(X0, X1, max_cells=200000):
    """
    Walk all attainable cells under a full PD metric.
    """
    seed = seed_cell(X0, X1)
    visited, attain = set(), []
    q = deque([seed])
    while q:
        if len(visited) >= max_cells:
            raise RuntimeError("interior cell walk exceeded max_cells = %d admissible cells; "
                               "the result would not be the exact optimum" % max_cells)
        A = q.popleft()
        if A in visited:
            continue
        visited.add(A)
        ok, _ = admissible(list(A), X0, X1)
        if not ok or not pd_attainable(A, X0):
            continue
        attain.append(A)
        for j in range(X0.shape[1]):
            if j in A:
                continue
            A2 = flip_neighbor(A, j, X0, X1)
            if A2 is not None and A2 not in visited:
                q.append(A2)
    return attain


def solve_interior_pd_exhaustive(X0, X1, Z1, Z0):
    """
    Evaluate every cell found by the PD-attainable cell walk.
    """
    cells = walk_attainable_cells_pd(X0, X1)
    best = None
    for A in cells:
        w = barycentric(list(A), X0, X1)
        wfull = np.zeros(X0.shape[1])
        wfull[list(A)] = w
        loss = float(np.sum((Z1 - Z0 @ wfull) ** 2))
        if best is None or loss < best["loss"]:
            best = dict(weights=wfull, loss=loss, cell=A, n_attainable=len(cells))
    return best


def solve_interior_pd(X0, X1, Z1, Z0, max_cells=200000):
    """
    Walk admissible cells and test PD attainability when a cell improves the loss.
    """
    J = X0.shape[1]

    from .interior_cellwalk import solve_interior_cellwalk
    diag = solve_interior_cellwalk(X0, X1, Z1, Z0)
    best = dict(weights=diag["weights"], loss=diag["loss"], cell=diag["cell"])
    U = diag["loss"]

    seed = seed_cell(X0, X1)
    visited = set()
    q = deque([seed])
    n_walked = n_sdp = 0
    while q:
        if len(visited) >= max_cells:
            raise RuntimeError("interior cell walk exceeded max_cells = %d admissible cells; "
                               "the result would not be the exact optimum" % max_cells)
        A = q.popleft()
        if A in visited:
            continue
        visited.add(A)
        luf = cell_factor(list(A), X0)              # one factorization, reused below
        ok, wA = admissible(list(A), X0, X1, luf)
        if not ok:
            continue
        n_walked += 1
        for A2 in cell_neighbors(A, X0, X1, luf):
            if A2 not in visited:
                q.append(A2)
        wf = np.zeros(J)
        wf[list(A)] = wA
        loss = float(np.sum((Z1 - Z0 @ wf) ** 2))
        if loss >= U - 1e-12:
            continue
        if theta_diag(list(A), X0, luf) > TOL:      # diagonal => PD-attainable, no SDP
            attainable = True
        else:
            n_sdp += 1
            attainable = theta_pd(list(A), X0, luf) > TOL
        if attainable:
            U = loss
            best = dict(weights=wf, loss=loss, cell=A)
    best["n_cells_walked"] = n_walked
    best["n_sdp_calls"] = n_sdp
    # Only improving cells are tested, so the total attainable count is unknown.
    best["n_attainable"] = None
    return best
