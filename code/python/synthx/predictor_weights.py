"""
Compute predictor weights that rationalize a given donor-weight vector.
"""

import itertools
import warnings
from collections import Counter

import cvxpy as cp
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

TOL = 1e-9


def _donor_geometry(X1, X0, w):
    """
    Residual and donor partition.
    """
    J = X0.shape[1]
    r = X1 - X0 @ w
    active = [j for j in range(J) if w[j] > 1e-7]
    inactive = [j for j in range(J) if w[j] <= 1e-7]
    if not active:
        raise ValueError("Donor weights must contain a positive weight.")
    return r, active, inactive


def _build_pd_optimality_conditions(X0, r, active, inactive):
    """
    Build V, the trace and active-donor equalities, and inactive-donor products.

    Each product is (X_j - X_j0)' V r and must be nonpositive.
    """
    K = X0.shape[0]
    j0 = active[0]
    V = cp.Variable((K, K), PSD=True)
    eqs = [cp.trace(V) == 1]
    for j in active[1:]:
        d = X0[:, j] - X0[:, j0]
        eqs.append(d @ V @ r == 0)
    products = ((j, (X0[:, j] - X0[:, j0]) @ V @ r) for j in inactive)
    return V, eqs, products


def _reported_number(value, decimals):
    return float(value) if decimals is None else round(float(value), decimals)


def _check_predictor_solution(problem):
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError("Predictor-weight optimization failed: " + problem.status)
    if problem.value is None or not np.isfinite(problem.value):
        raise RuntimeError("Predictor-weight optimization returned no finite optimum.")


def compute_pd_logdet_and_ranges(X1, X0, w, names, *, decimals=4, ranges=True):
    """
    Compute the PD log-determinant maximizer and predictor-weight ranges.

    Use decimals=None for full precision and ranges=False to skip the ranges.
    """
    K = X0.shape[0]
    r, active, inactive = _donor_geometry(X1, X0, w)
    V, eqs, products = _build_pd_optimality_conditions(X0, r, active, inactive)
    cons = eqs + [product <= 0 for _, product in products]

    prob = cp.Problem(cp.Maximize(cp.log_det(V)), cons)
    prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-12,
               tol_gap_rel=1e-12, tol_feas=1e-12)
    _check_predictor_solution(prob)
    Vstar = V.value.copy()
    eigs = np.linalg.eigvalsh(Vstar)
    out = dict(status=prob.status,
               logdet=float(prob.value),
               min_eig=float(eigs.min()),
               diag={names[k]: _reported_number(Vstar[k, k], decimals)
                     for k in range(K)})
    if not ranges:
        return Vstar, out

    identified_ranges = {}
    inaccurate = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # cvxpy's 'inaccurate' warning is reported below instead
        for k in range(K):
            lo = cp.Problem(cp.Minimize(V[k, k]), cons)
            lo.solve(solver=cp.CLARABEL)
            _check_predictor_solution(lo)
            hi = cp.Problem(cp.Maximize(V[k, k]), cons)
            hi.solve(solver=cp.CLARABEL)
            _check_predictor_solution(hi)
            if "inaccurate" in (lo.status + hi.status):
                inaccurate.append(names[k])
            lower = float(lo.value)
            if decimals is not None and abs(lower) < 0.5 * 10 ** (-decimals):
                lower = abs(lower)  # Avoid displaying negative zero.
            identified_ranges[names[k]] = [
                _reported_number(lower, decimals),
                _reported_number(hi.value, decimals),
            ]
    if inaccurate:
        warnings.warn(
            "CLARABEL reported an inaccurate solution for the identified range of: "
            + ", ".join(map(str, inaccurate)),
            RuntimeWarning,
            stacklevel=2,
        )
    out["ranges"] = identified_ranges
    return Vstar, out


def analytic_center(X1, X0, w, names, *, decimals=4):
    """
    Maximize log det V plus the log of each slack that can be positive.

    Slacks that are always zero are excluded. Use decimals=None for full precision.
    """
    K = X0.shape[0]
    r, active, inactive = _donor_geometry(X1, X0, w)
    j0 = active[0]
    V, eqs, products = _build_pd_optimality_conditions(X0, r, active, inactive)
    slack = {j: -product for j, product in products}
    ineqs = [slack[j] >= 0 for j in inactive]

    # Check whether all inactive slacks can be positive together.
    strict = []
    if inactive:
        t = cp.Variable()
        joint = cp.Problem(cp.Maximize(t),
                           eqs + [slack[j] >= t for j in inactive])
        joint.solve(solver=cp.CLARABEL, tol_gap_abs=1e-12,
                    tol_gap_rel=1e-12, tol_feas=1e-12)
        _check_predictor_solution(joint)
        if joint.value > TOL:
            strict = list(inactive)
        else:
            # A slack that is always zero cannot enter the log barrier.
            for j in inactive:
                pj = cp.Problem(cp.Maximize(slack[j]), eqs + ineqs)
                pj.solve(solver=cp.CLARABEL, tol_gap_abs=1e-12,
                         tol_gap_rel=1e-12, tol_feas=1e-12)
                _check_predictor_solution(pj)
                if pj.value > TOL:
                    strict.append(j)

    obj = cp.log_det(V)
    if strict:
        obj = obj + cp.sum(cp.hstack([cp.log(slack[j]) for j in strict]))
    prob = cp.Problem(cp.Maximize(obj), eqs + ineqs)
    prob.solve(solver=cp.CLARABEL)
    _check_predictor_solution(prob)
    Vc = V.value
    eigs = np.linalg.eigvalsh(Vc)
    slacks = {j: float(-((X0[:, j] - X0[:, j0]) @ Vc @ r))
              for j in inactive}
    return Vc, dict(
        status=prob.status,
        n_inactive=len(inactive),
        n_implicit_equalities=len(inactive) - len(strict),
        min_eig=float(eigs.min()),
        min_strict_slack=(min(slacks[j] for j in strict)
                          if strict else None),
        diag={names[k]: _reported_number(Vc[k, k], decimals) for k in range(K)},
    )


def diag_sets(X0, X1, w):
    r, active, inactive = _donor_geometry(X1, X0, w)
    j0 = active[0]
    A_eq = (X0[:, active[1:]] - X0[:, [j0]]).T * r
    S = (X0[:, [j0]] - X0[:, inactive]).T * r
    return r, inactive, A_eq, S


def lp_max(c, A_eq, S, K):
    """
    Maximize c'd over the diagonal predictor weights that rationalize w.
    """
    Aeq = np.vstack([A_eq, np.ones((1, K))])
    beq = np.zeros(Aeq.shape[0])
    beq[-1] = 1.0
    res = linprog(-c, A_ub=-S, b_ub=np.zeros(S.shape[0]), A_eq=Aeq, b_eq=beq,
                  bounds=(0, None), method="highs")
    if res.status != 0:
        raise RuntimeError("Diagonal predictor-weight LP failed: " + res.message)
    return -float(res.fun)


def diag_centers(X0, X1, w, names, *, decimals=4):
    """
    Log-det maximizer and analytic center of the diagonal rationalizers.

    Exclude weights and slacks that must be zero from the log barriers.
    Use decimals=None for full precision.
    """
    K = X0.shape[0]
    r, inactive, A_eq, S = diag_sets(X0, X1, w)
    support = [k for k in range(K) if lp_max(np.eye(K)[k], A_eq, S, K) > TOL]
    forced = [k for k in range(K) if k not in support]
    smax = [lp_max(S[m], A_eq, S, K) for m in range(len(inactive))]
    strict = [m for m in range(len(inactive)) if smax[m] > TOL]
    implicit = [m for m in range(len(inactive)) if smax[m] <= TOL]

    d = cp.Variable(K, nonneg=True)
    cons = [cp.sum(d) == 1]
    if len(A_eq):
        cons.append(A_eq @ d == 0)
    if len(S):
        cons.append(S @ d >= 0)
    cons += [d[k] == 0 for k in forced] + [S[m] @ d == 0 for m in implicit]
    out = {}
    objectives = {"logdet": cp.sum(cp.log(d[support])),
                  "analytic_center": cp.sum(cp.log(d[support])) + (cp.sum(cp.log(S[strict] @ d)) if strict else 0)}
    for label, obj in objectives.items():
        prob = cp.Problem(cp.Maximize(obj), cons)
        prob.solve(solver=cp.CLARABEL)
        _check_predictor_solution(prob)
        dv = np.asarray(d.value).ravel()
        out[label] = {names[k]: _reported_number(dv[k], decimals) for k in range(K)}
        out[label + "_status"] = prob.status
    out["support"] = [names[k] for k in support]
    out["forced_zero"] = [names[k] for k in forced]
    out["n_inactive_donors"] = len(inactive)
    out["n_implicit_equalities"] = len(implicit)
    out["residual_on_support"] = {names[k]: float(r[k]) for k in support}
    return out


def standardize(X1, X0):
    """
    Divide each predictor by its sample standard deviation across all units.
    """
    X1 = np.asarray(X1, float).reshape(-1)
    X0 = np.asarray(X0, float)
    sd = np.column_stack([X1, X0]).std(axis=1, ddof=1)
    sd = np.where(sd > 0, sd, 1.0)
    return X1 / sd, X0 / sd[:, None], sd


def diag_rationalizers(X1, X0, w, standardized=False):
    """
    Return full-precision diagonal centers, ranges, and retained predictors.
    """
    X1 = np.asarray(X1, float).reshape(-1)
    X0 = np.asarray(X0, float)
    if not standardized:
        X1, X0, _ = standardize(X1, X0)
    K = X0.shape[0]
    names = list(range(K))
    centers = diag_centers(X0, X1, w, names, decimals=None)
    _, _, A_eq, S = diag_sets(X0, X1, w)
    lo = np.array([-lp_max(-row, A_eq, S, K) for row in np.eye(K)])
    hi = np.array([lp_max(row, A_eq, S, K) for row in np.eye(K)])
    return dict(
        logdet=np.array([centers["logdet"][k] for k in names]),
        center=np.array([centers["analytic_center"][k] for k in names]),
        lo=lo,
        hi=hi,
        support=np.array(centers["support"], dtype=int),
        implicit_eq=centers["n_implicit_equalities"],
    )


def pd_rationalizers(X1, X0, w, standardized=False, ranges=True):
    """
    Return full-precision PD centers and optional diagonal-entry ranges.
    """
    X1 = np.asarray(X1, float).reshape(-1)
    X0 = np.asarray(X0, float)
    if not standardized:
        X1, X0, _ = standardize(X1, X0)
    names = list(range(X0.shape[0]))
    Vstar, logdet = compute_pd_logdet_and_ranges(
        X1, X0, w, names, decimals=None, ranges=ranges,
    )
    center, center_info = analytic_center(X1, X0, w, names, decimals=None)
    result = dict(
        Vstar=Vstar,
        center=center,
        logdet_diag=np.diag(Vstar).copy(),
        center_diag=np.diag(center).copy(),
        min_eig=logdet["min_eig"],
        center_min_eig=center_info["min_eig"],
        implicit_eq=center_info["n_implicit_equalities"],
    )
    if ranges:
        result["lo"] = np.array([logdet["ranges"][k][0] for k in names])
        result["hi"] = np.array([logdet["ranges"][k][1] for k in names])
    return result


def qhull_facets(X0, tol=1e-9):
    """
    Find the donors on each distinct supporting plane of the hull.
    """
    hull = ConvexHull(X0.T)
    eqs = {}
    for eq in hull.equations:
        eqs[tuple(np.round(eq, 9))] = np.array(eq)
    out = []
    pts = X0.T
    for eq in eqs.values():
        a, b = eq[:-1], eq[-1]
        on = np.where(np.abs(pts @ a + b) <= tol * max(1.0, abs(b)))[0]
        out.append((frozenset(int(i) for i in on), a, b))
    return out, hull


def faces_below(facet_sets):
    faces = set()
    for fs in facet_sets:
        verts = sorted(fs)
        for k in range(1, len(verts) + 1):
            for sub in itertools.combinations(verts, k):
                faces.add(frozenset(sub))
    return faces


def counts_by_dim(faces, K):
    c = Counter(len(f) - 1 for f in faces)
    return [int(c.get(d, 0)) for d in range(K)]


def count_faces(X1, X0, standardized=False, progress=None):
    """
    Count attainable, visible, and all faces in dimensions 0 through K - 1.

    As in the paper, faces must have affinely independent donors.
    """
    from .visible_facet_solver import feasible_pieces, visible_facets

    X1 = np.asarray(X1, float).reshape(-1)
    X0 = np.asarray(X0, float)
    if not standardized:
        X1, X0, _ = standardize(X1, X0)
    K = X0.shape[0]
    vf = visible_facets(X0, X1)
    if K == 1:
        # A one-dimensional hull has its minimum and maximum as facets.
        facet_sets = [frozenset([int(np.argmin(X0[0]))]),
                      frozenset([int(np.argmax(X0[0]))])]
        vis_q = set(vf)
    else:
        facets, _ = qhull_facets(X0)
        facet_sets = [face for face, _, _ in facets]
        vis_q = {face for face, a, b in facets if a @ X1 + b > TOL}
    all_faces = faces_below(facet_sets)
    visible_faces = faces_below(vis_q)
    attainable = set()
    for i, face in enumerate(sorted(visible_faces, key=len, reverse=True)):
        normal = vf[face][0] if face in vf else None
        if feasible_pieces(face, X0, X1, normal):
            attainable.add(face)
        if progress is not None and (i + 1) % 250 == 0:
            progress(i + 1, len(visible_faces))
    return np.column_stack([
        counts_by_dim(attainable, K),
        counts_by_dim(visible_faces, K),
        counts_by_dim(all_faces, K),
    ])
