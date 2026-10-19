"""
Find attainable interior cells by walking between neighboring cells.
"""

from __future__ import annotations

from collections import deque
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import lu_solve

from .interior_vsearch import admissible, theta_diag, barycentric, cell_factor

TOL = 1e-9
LP = dict(method="highs")


def affine_dependence(B, X0):
    """
    Return the affine dependence of K+2 donors, aligned with B.
    """
    M = np.vstack([X0[:, B], np.ones(len(B))])          # (K+1) x (K+2)
    _, _, Vt = np.linalg.svd(M)
    return Vt[-1]                                         # right null vector


def _det_perturb(J, attempt):
    """
    Return a deterministic direction for breaking ties in the seed LP.
    """
    phi = (1.0 + 5.0 ** 0.5) / 2.0
    psi = 2.0 ** 0.5
    j = np.arange(1, J + 1)
    return np.mod(j * phi + attempt * psi, 1.0) - 0.5


def seed_cell(X0, X1):
    """
    Find a full-dimensional cell containing X1 using the lifted-distance LP.
    """
    K, J = X0.shape
    A_eq = np.vstack([X0, np.ones(J)])
    b_eq = np.append(X1, 1.0)
    base = np.sum(X0 ** 2, axis=0)
    for attempt in range(6):
        c = base + (0.0 if attempt == 0 else 1e-6 * _det_perturb(J, attempt))
        res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), **LP)
        if res.status != 0:
            continue
        supp = tuple(sorted(j for j in range(J) if res.x[j] > 1e-7))
        if len(supp) == K + 1:
            return supp
    raise RuntimeError("could not seed a full-dimensional interior cell")


def flip_neighbor(A, j, X0, X1):
    """
    Replace one cell vertex with donor j, keeping a cell that contains X1.
    """
    B = list(A) + [j]
    lam = affine_dependence(B, X0)
    sj = lam[-1]                                          # coefficient of j
    if abs(sj) < TOL:
        return None
    for idx in range(len(B)):
        if lam[idx] * sj < -TOL:                          # opposite sign to j
            C = [B[t] for t in range(len(B)) if t != idx]
            try:
                w = barycentric(C, X0, X1)
            except np.linalg.LinAlgError:
                continue
            if (w >= -1e-9).all():
                return tuple(sorted(C))
    return None


def cell_neighbors(A, X0, X1, luf):
    """
    Find neighboring cells using the stored factorization of cell A.
    """
    A = list(A)
    nA = len(A)
    nbrs = []
    for j in range(X0.shape[1]):
        if j in A:
            continue
        beta = lu_solve(luf, np.append(X0[:, j], 1.0))    # barycentric of j in A
        bb = float(beta @ beta) + 1.0                     # ||(beta,-1)||^2
        if bb > 1.0 / TOL ** 2:                           # abs(sj) < TOL: no circuit
            continue
        for idx in range(nA):
            if beta[idx] > TOL * bb:                      # == lam[idx]*sj < -TOL
                C = [A[t] for t in range(nA) if t != idx] + [j]
                try:
                    w = barycentric(C, X0, X1)
                except np.linalg.LinAlgError:
                    continue
                if (w >= -1e-9).all():
                    nbrs.append(tuple(sorted(C)))
                    break
    return nbrs


def walk_attainable_cells(X0, X1, max_cells=200000):
    """
    Walk the diagonal-attainable cells, factoring each cell once.
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
        luf = cell_factor(list(A), X0)
        ok, _ = admissible(list(A), X0, X1, luf)
        if not ok or theta_diag(list(A), X0, luf) <= TOL:
            continue
        attain.append(A)
        for A2 in cell_neighbors(A, X0, X1, luf):
            if A2 not in visited:
                q.append(A2)
    return attain


def solve_interior_cellwalk(X0, X1, Z1, Z0):
    """
    Select the attainable interior cell with the smallest validation loss.
    """
    cells = walk_attainable_cells(X0, X1)
    best = None
    for A in cells:
        w = barycentric(list(A), X0, X1)
        wfull = np.zeros(X0.shape[1])
        wfull[list(A)] = w
        loss = float(np.sum((Z1 - Z0 @ wfull) ** 2))
        if best is None or loss < best["loss"]:
            best = dict(weights=wfull, loss=loss, cell=A,
                        n_attainable=len(cells))
    if best is None:
        raise RuntimeError("no attainable interior cell found for the diagonal class")
    return best
