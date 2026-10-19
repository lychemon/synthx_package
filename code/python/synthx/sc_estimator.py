"""
Choose the exterior, interior, or fixed-metric donor-weight solver.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import linprog

from .visible_facet_solver import solve_exterior_hullfree as solve_exterior
from .visible_facet_solver import solve_exterior_pd
from .visible_facet_solver import _affine_rank
from .interior_cellwalk import solve_interior_cellwalk
from .interior_pd import solve_interior_pd
from .simplex_solver import qp_simplex
from .gpcheck import general_position_check

LP_OPTS = dict(method="highs")


def inverse_variance_V(X0, X1):
    """
    Return the inverse-variance diagonal metric on pooled predictors.
    """
    pooled = np.column_stack([X0, X1.reshape(-1, 1)])   # K x (J+1)
    v = np.var(pooled, axis=1)
    v[v <= 0] = 1.0
    return np.diag(1.0 / v)


def _as_pd_V(V, K):
    """
    Validate a positive diagonal vector or a symmetric positive definite matrix.
    """
    V = np.asarray(V, float)
    if not np.isfinite(V).all():
        raise ValueError("V must contain finite values")
    if V.ndim == 1:
        if V.shape[0] != K:
            raise ValueError(
                f"diagonal V has length {V.shape[0]}, expected K={K}")
        if np.any(V <= 0):
            raise ValueError(
                "diagonal V (a vector) must be positive definite: every "
                "entry must be > 0")
        return np.diag(V)
    if V.ndim == 2:
        if V.shape != (K, K):
            raise ValueError(
                f"V has shape {V.shape}, expected a {(K, K)} matrix or a "
                f"length-{K} vector")
        if not np.allclose(V, V.T, atol=1e-10):
            raise ValueError("V must be symmetric")
        Vs = 0.5 * (V + V.T)
        if np.linalg.eigvalsh(Vs).min() <= 0:
            raise ValueError(
                "V must be positive definite (all eigenvalues > 0)")
        return Vs
    raise ValueError(
        f"V must be a length-{K} vector or a {K} x {K} matrix")


def solve_fixed_V(X0, X1, Z1, Z0, V):
    """
    Minimize predictor distance under V, then report validation loss.
    """
    K = X0.shape[0]
    Vm = _as_pd_V(V, K)
    L = np.linalg.cholesky(Vm)
    # Remove common offsets and the metric's arbitrary overall scale.
    differences = L.T @ (X0 - X1[:, None])
    scale = float(np.abs(differences).max())
    if scale > 0:
        differences = differences / scale
    _, w = qp_simplex(np.zeros(K), differences)
    w = np.clip(np.asarray(w, float), 0, None)
    w[w < 1e-7] = 0.0
    w = w / w.sum()
    r = Z1 - Z0 @ w
    return float(r @ r), w, dict(V=Vm)


def in_hull(X1, X0, tol=1e-9):
    """
    Test whether X1 is strictly inside the donor hull; also return the LP margin.
    """
    K, J = X0.shape
    c = np.zeros(J + 1)
    c[-1] = -1.0                       # maximize t
    A_eq = np.zeros((K + 1, J + 1))
    A_eq[:K, :J] = X0
    A_eq[K, :J] = 1.0
    b_eq = np.concatenate([X1, [1.0]])
    A_ub = np.zeros((J, J + 1))        # -w + t <= 0  ->  w >= t
    A_ub[:, :J] = -np.eye(J)
    A_ub[:, J] = 1.0
    b_ub = np.zeros(J)
    bounds = [(0, None)] * J + [(None, None)]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds, **LP_OPTS)
    if res.status == 2:
        return False, 0.0
    if res.status != 0:
        raise RuntimeError("Hull-membership LP failed: " + res.message)
    t = -res.fun
    return t > tol, float(t)


def _donor_conditioning(X0):
    """
    Return the condition number after centering and scaling donor predictors.
    """
    C = X0 - X0.mean(axis=1, keepdims=True)
    sd = C.std(axis=1)
    sd = np.where(sd > 0, sd, 1.0)
    sv = np.linalg.svd(C / sd[:, None], compute_uv=False)
    if sv.size == 0 or sv[-1] <= 0:
        return np.inf
    return float(sv[0] / sv[-1])


def _baseline_guard(kappa, X0, X1, Z1, Z0, loss, w, info, V_class, slack, regime):
    """
    For poorly conditioned data, compare the fit with inverse-variance weights.
    """
    loss = float(loss)
    if kappa > 1e5:
        try:
            loss_iv, w_iv, _ = solve_fixed_V(X0, X1, Z1, Z0,
                                             inverse_variance_V(X0, X1))
        except Exception:
            loss_iv = None
        if loss_iv is not None and loss > loss_iv + 1e-6 * max(1.0, abs(loss_iv)):
            warnings.warn(
                "The %s solve gave loss %.6g; the inverse-variance baseline gave "
                "%.6g. Using baseline weights because the predictor geometry "
                "is poorly conditioned." % (V_class, loss, loss_iv),
                stacklevel=3)
            w, loss, info = w_iv, float(loss_iv), dict(info, baseline_fallback=True)
            info.pop("general_position", None)
    return dict(weights=w, loss=loss, regime=regime, V_class=V_class,
                interior_slack=slack, info=info)


def _solve_pd_lowdim(X0, X1, Z1, Z0, slack, tol=1e-7):
    """
    Solve the PD problem when donors span fewer than K dimensions.
    """
    c0 = X0.mean(axis=1)
    C = X0 - c0[:, None]
    U, sv, _ = np.linalg.svd(C, full_matrices=True)
    r = int((sv > tol * sv[0]).sum()) if sv.size and sv[0] > 0 else 0
    B = U[:, :r]
    T = B.T @ C
    tau = B.T @ (X1 - c0)
    offspan = float(np.linalg.norm((X1 - c0) - B @ tau))
    outside_span = offspan > tol * max(1.0, np.linalg.norm(X1 - c0))
    if outside_span or r == 0:
        # Outside the donor span, a PD metric can rationalize any donor mixture.
        # Identical donors also leave every mixture available.
        _, w = qp_simplex(Z1, Z0)
        residual = Z1 - Z0 @ w
        return dict(weights=w, loss=float(residual @ residual),
                    regime="exterior" if outside_span else "interior",
                    V_class="pd", interior_slack=slack,
                    info=dict(reduced_rank=r, offspan_norm=offspan,
                              path="pd-offspan" if outside_span else "pd-identical"))
    warnings.warn(
        "Solving the PD problem in the donors' %d-dimensional affine span." % r,
        stacklevel=2)
    inside, _ = in_hull(tau, T)
    if inside:
        sol = solve_interior_pd(T, tau, Z1, Z0)
        w, loss = sol["weights"], sol["loss"]
        info = {k: v for k, v in sol.items() if k not in ("weights", "loss")}
        regime = "interior"
    else:
        loss, w, info = solve_exterior_pd(T, tau, Z1, Z0)
        regime = "exterior"
    return dict(weights=w, loss=float(loss), regime=regime, V_class="pd",
                interior_slack=slack,
                info=dict(info, reduced_rank=r, offspan_norm=offspan,
                          path="pd-reduced"))


def estimate(X1, X0, Z1, Z0, V_class="diagonal", V=None, donor_labels=None):
    """
    Estimate donor weights from predictors X and validation outcomes Z.

    X1 and Z1 describe the treated unit; columns of X0 and Z0 describe donors.
    V_class selects diagonal or full PD optimization. Supplying V fixes the metric.
    Return weights, loss, regime, V_class, interior_slack, and solver information.
    """
    X1 = np.asarray(X1, float)
    X0 = np.asarray(X0, float)
    Z1 = np.asarray(Z1, float)
    Z0 = np.asarray(Z0, float)
    if (X0.ndim != 2 or Z0.ndim != 2 or X1.shape != (X0.shape[0],)
            or Z1.shape != (Z0.shape[0],) or X0.shape[1] != Z0.shape[1]
            or min(X0.shape) == 0 or Z0.shape[0] == 0):
        raise ValueError("Use treated vectors and donor matrices with matching rows and donor columns")
    if not all(np.isfinite(a).all() for a in (X1, X0, Z1, Z0)):
        raise ValueError("Predictors and validation outcomes must contain finite values")

    if V is not None:
        loss, w, info = solve_fixed_V(X0, X1, Z1, Z0, V)
        interior, slack = in_hull(X1, X0)
        return dict(weights=w, loss=loss,
                    regime="interior" if interior else "exterior",
                    V_class="fixed", interior_slack=slack, info=info)

    if V_class not in ("diagonal", "pd"):
        raise ValueError(f"V_class must be 'diagonal' or 'pd', got {V_class!r}")

    n_dropped = 0
    if V_class == "diagonal":
        spread = np.ptp(X0, axis=1)
        scale = np.maximum(np.abs(X0).max(axis=1), 1.0)
        keep = spread > 1e-12 * scale
        n_dropped = int((~keep).sum())
        if n_dropped:
            X0, X1 = X0[keep], X1[keep]
            warnings.warn("Dropped %d constant donor predictors for diagonal V." % n_dropped,
                          stacklevel=2)
    if X0.shape[0] == 0:
        raise ValueError(
            "All donor predictors are constant. Supply a fixed metric with V=.")

    # Use the same predictor scale as the paper. A supplied V keeps its original scale.
    sd = np.column_stack([X1, X0]).std(axis=1, ddof=1)
    sd[sd == 0] = 1.0
    X1, X0 = X1 / sd, X0 / sd[:, None]
    interior, slack = in_hull(X1, X0)
    if _affine_rank(X0) < X0.shape[0]:
        if V_class == "pd":
            result = _solve_pd_lowdim(X0, X1, Z1, Z0, slack)
            result["info"]["n_dropped_constant"] = n_dropped
            return result
        raise ValueError(
            "Donor predictors are collinear after removing constant rows. "
            "The diagonal solver requires a full-dimensional hull; supply a fixed V "
            "or revise the predictors. Removing a collinear row changes the diagonal problem.")

    kappa = _donor_conditioning(X0)
    if kappa > 1e5:
        warnings.warn(
            "Donor predictors are nearly collinear (condition number %.1e); "
            "the result may be imprecise." % kappa,
            stacklevel=2)

    if not interior:
        solver = solve_exterior if V_class == "diagonal" else solve_exterior_pd
        loss, w, info = solver(X0, X1, Z1, Z0)
        face = info.get("face")
        if face is None:
            face = tuple(np.flatnonzero(w > 1e-7))
        gp = general_position_check(X0, X1, face, info.get("visible_facets"),
                                    diagonal=(V_class == "diagonal"), labels=donor_labels)
        info = {k: v for k, v in info.items() if k != "visible_facets"}
        return _baseline_guard(kappa, X0, X1, Z1, Z0, loss, w,
                               dict(info, n_dropped_constant=n_dropped, general_position=gp),
                               V_class, slack, "exterior")

    sol = (solve_interior_cellwalk if V_class == "diagonal"
           else solve_interior_pd)(X0, X1, Z1, Z0)
    info = {k: v for k, v in sol.items() if k not in ("weights", "loss")}
    info["n_dropped_constant"] = n_dropped
    return _baseline_guard(kappa, X0, X1, Z1, Z0, sol["loss"], sol["weights"],
                           info, V_class, slack, "interior")


def predictor_balance(X1, X0, w):
    """
    Return the treated and synthetic predictor vectors.
    """
    return X1, X0 @ np.asarray(w, float)
