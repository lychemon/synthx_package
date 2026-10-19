"""
Fit outcomes using nonnegative donor weights that sum to one.
"""

from __future__ import annotations

import cvxpy as cp
import numpy as np


def affine_bound(Z1, Zc):
    """
    Use the affine fit to bound the simplex loss.
    """
    target = Z1 - Zc[:, 0]
    directions = Zc[:, 1:] - Zc[:, :1]
    if directions.shape[1] == 0:
        weights = np.ones(1)
    else:
        coefficients, _, rank, _ = np.linalg.lstsq(directions, target, rcond=1e-10)
        if rank < directions.shape[1]:
            # Use zero when the donor directions are numerically dependent.
            return 0.0
        weights = np.concatenate(([1.0 - coefficients.sum()], coefficients))
    return simplex_lower_bound(Z1, Zc, weights)


def simplex_lower_bound(Z1, Zc, weights):
    """
    Bound the simplex optimum using the loss and gradient at any weights.
    """
    residual = Z1 - Zc @ weights
    loss = float(residual @ residual)
    gradient = -2.0 * Zc.T @ residual
    # The tangent plane lies below the loss; its simplex minimum is at a donor.
    bound = loss + float(gradient.min() - gradient @ weights)
    roundoff = 1e-10 * max(1.0, loss, float(np.abs(gradient).max()),
                           abs(float(gradient @ weights)))
    return max(0.0, bound - roundoff)


def qp_simplex(Z1, Zc, G=None, h=None, tol=1e-12, maxiter=200,
               A_extra=None, b_extra=None):
    """
    Minimize ||Z1 - Zc @ b||^2 with nonnegative weights that sum to one.
    Optionally require G @ b <= h and A_extra @ b = b_extra.

    Return (loss, weights), or (None, None) if no feasible weights exist.
    Raise an error if Clarabel does not converge or returns infeasible weights.
    """
    m = Zc.shape[1]
    b = cp.Variable(m)
    constraints = [b >= 0, cp.sum(b) == 1]
    inequalities = G is not None and len(G) > 0
    equalities = A_extra is not None and len(A_extra) > 0
    if inequalities:
        G, h = np.asarray(G, float), np.asarray(h, float)
        scale_G = np.maximum(1.0, np.maximum(np.abs(G).max(axis=1), np.abs(h)))
        G, h = G / scale_G[:, None], h / scale_G
        constraints.append(G @ b <= h)
    if equalities:
        A_extra, b_extra = np.asarray(A_extra, float), np.asarray(b_extra, float)
        scale_A = np.maximum(1.0, np.maximum(np.abs(A_extra).max(axis=1),
                                            np.abs(b_extra)))
        A_extra, b_extra = A_extra / scale_A[:, None], b_extra / scale_A
        constraints.append(A_extra @ b == b_extra)

    # Center and scale the outcomes without changing the minimizing weights.
    differences = Zc - Z1[:, None]
    scale = max(1.0, float(np.abs(differences).max()))
    outcome_rows = differences / scale
    if outcome_rows.shape[0] > m:
        # QR removes extra outcome rows while preserving the squared loss.
        _, outcome_rows = np.linalg.qr(outcome_rows, mode="reduced")
    problem = cp.Problem(cp.Minimize(cp.sum_squares(outcome_rows @ b)),
                         constraints)
    try:
        problem.solve(solver=cp.CLARABEL, max_iter=maxiter,
                      tol_gap_abs=tol, tol_gap_rel=tol, tol_feas=tol)
    except cp.error.SolverError as error:
        raise RuntimeError("The donor-weight QP solver failed.") from error
    if problem.status == cp.INFEASIBLE:
        return None, None
    if problem.status != cp.OPTIMAL:
        raise RuntimeError("The donor-weight QP did not converge: " + problem.status)

    weights = np.asarray(b.value).reshape(-1)
    violations = [abs(weights.sum() - 1.0), -weights.min()]
    if inequalities:
        violations.append(float(np.max(G @ weights - h)))
    if equalities:
        violations.append(float(np.max(np.abs(A_extra @ weights - b_extra))))
    if not np.all(np.isfinite(weights)) or max(violations) > max(1e-9, 10 * tol):
        raise RuntimeError("The donor-weight QP returned infeasible weights.")
    residual = Z1 - Zc @ weights
    return float(residual @ residual), weights
