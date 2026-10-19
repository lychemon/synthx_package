"""
Find visible facets and solve diagonal and positive-definite donor-weight problems.
"""

from __future__ import annotations

import heapq
import itertools
import math
from collections import deque

import numpy as np
from scipy.optimize import linprog

from .simplex_solver import qp_simplex, affine_bound, simplex_lower_bound

LP = dict(method="highs")
TOL = 1e-9


class LowerDimHull(RuntimeError):
    """
    The donor hull is (numerically) lower-dimensional than R^K, so no
    full-dimensional visible facet exists.
    """


# ----------------------------------------------------------------------
# Facet geometry (linear algebra only)
# ----------------------------------------------------------------------
def interior_point(X0):
    return X0.mean(axis=1)


def facet_normal(X0, S, interior):
    """
    Outward unit normal n and offset b of the hyperplane through donors S
    (|S| = K), oriented so the hull interior satisfies n.x < b. None if the
    points are affinely degenerate.
    """
    pts = X0[:, S]
    base = pts[:, 0]
    E = pts[:, 1:] - base[:, None]          # K x (K-1)
    U, sv, _ = np.linalg.svd(E)
    if sv.size > 0 and sv[-1] < 1e-12:
        return None
    n = U[:, -1]                            # orthogonal to span(E)
    b = n @ base
    if n @ interior > b:
        n, b = -n, -b
    return n, b


def all_donors_in_halfspace(n, b, X0):
    """
    Return whether every donor satisfies n @ x <= b within tolerance 1e-7.
    """
    return float((n @ X0).max()) <= b + 1e-7


def visible(n, b, X1):
    return n @ X1 > b + TOL


# ----------------------------------------------------------------------
# Find a starting facet visible from X1.
# First try where the ray from the hull interior toward X1 leaves the hull.
# If it leaves at a vertex or ridge, find a facet containing that point.
# If needed, try a separating hyperplane to find candidate facet vertices.
# Raise LowerDimHull if the donor hull has dimension less than K.
# ----------------------------------------------------------------------
def _affine_rank(X0, tol=1e-7):
    """
    Dimension of the affine hull of the donors (rank of the centered set).
    """
    C = X0 - X0.mean(axis=1, keepdims=True)
    if C.shape[1] == 0:
        return 0
    sv = np.linalg.svd(C, compute_uv=False)
    return int((sv > tol * sv[0]).sum()) if sv[0] > 0 else 0


def _pick_facet_vertices(X0, cand, interior):
    """
    Choose K donors from `cand` that define a hyperplane.
    Return (indices, (normal, offset)), or None if unsuccessful.
    """
    K = X0.shape[0]
    cand = sorted(cand)
    if len(cand) < K:
        return None
    chosen = [cand[0]]
    base = X0[:, cand[0]]
    E = np.zeros((K, 0))
    for j in cand[1:]:
        E2 = np.hstack([E, (X0[:, j] - base)[:, None]])
        if np.linalg.matrix_rank(E2, tol=1e-7) > E.shape[1]:
            E = E2
            chosen.append(j)
        if len(chosen) == K:
            break
    if len(chosen) < K:
        return None
    fn = facet_normal(X0, sorted(chosen), interior)
    if fn is None:
        return None
    return frozenset(chosen), fn


def _separation_normal(X0, X1):
    """
    Outward normal (n, b) of a hyperplane separating X1 from conv(X0),
    via an LP:  max n.X1 - b  s.t. n.X0_j - b <= 0, -1 <= n <= 1.
    None if the LP fails.
    """
    K, J = X0.shape
    c = np.concatenate([-X1, [1.0]])
    A_ub = np.hstack([X0.T, -np.ones((J, 1))])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(J),
                  bounds=[(-1, 1)] * K + [(None, None)], **LP)
    if res.status != 0:
        return None
    return res.x[:K], res.x[K]


def _facet_at_point(X0, point, interior):
    """
    Find a facet containing a point on the hull boundary.
    """
    K, J = X0.shape
    A_eq = [point - interior]
    b_eq = [1.0]
    # Fix the normal's scale, then minimize its coordinates in order.
    # Keep earlier minima to resolve ties until the plane contains a facet.
    for k in range(K):
        objective = np.eye(K)[k]
        res = linprog(
            objective,
            A_ub=(X0 - point[:, None]).T, b_ub=np.zeros(J),
            A_eq=A_eq, b_eq=b_eq,
            bounds=(None, None), method="highs-ds",
        )
        if res.status != 0:
            return None
        n = res.x / np.linalg.norm(res.x)
        b = n @ point
        candidates = np.flatnonzero(np.abs(n @ X0 - b) <= 1e-7).tolist()
        facet = _pick_facet_vertices(X0, candidates, interior)
        if facet is not None:
            return facet
        A_eq.append(objective)
        b_eq.append(res.x[k])
    return None


def seed_visible_facet(X0, X1, interior):
    K, J = X0.shape
    if _affine_rank(X0) < K:                       # no full-dimensional facet
        raise LowerDimHull(
            "donor hull is lower-dimensional than R^%d (rank deficient); "
            "no full-dimensional visible facet exists" % K)

    # Find where the ray from the hull interior toward X1 leaves the hull.
    dirv = X1 - interior
    c = np.zeros(J + 1)
    c[-1] = -1.0
    A_eq = np.zeros((K + 1, J + 1))
    A_eq[:K, :J] = X0
    A_eq[:K, J] = -dirv
    A_eq[K, :J] = 1.0
    b_eq = np.concatenate([interior, [1.0]])
    res = linprog(c, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, None)] * J + [(0, None)], **LP)
    if res.status == 0:
        supp = [j for j in range(J) if res.x[j] > 1e-7]
        if len(supp) == K:                         # clean exit facet
            fn = facet_normal(X0, supp, interior)
            if fn and all_donors_in_halfspace(*fn, X0) and visible(*fn, X1):
                return frozenset(supp)
        if len(supp) >= K:                         # degenerate exit face: recover
            got = _pick_facet_vertices(X0, supp, interior)
            if (got is not None and all_donors_in_halfspace(*got[1], X0)
                    and visible(*got[1], X1)):
                return got[0]

        # The ray can leave through a vertex or ridge shared by several facets.
        exit_point = X0 @ res.x[:J]
        got = _facet_at_point(X0, exit_point, interior)
        if (got is not None and all_donors_in_halfspace(*got[1], X0)
                and visible(*got[1], X1)):
            return got[0]

    # Try the donors on a separating plane if the exit-point search failed.
    sep = _separation_normal(X0, X1)
    if sep is not None:
        n, b = sep
        cand = [j for j in range(J) if n @ X0[:, j] >= b - 1e-7]
        got = _pick_facet_vertices(X0, cand, interior)
        if (got is not None and all_donors_in_halfspace(*got[1], X0)
                and visible(*got[1], X1)):
            return got[0]

    raise RuntimeError("could not seed a visible facet")


# ----------------------------------------------------------------------
# Find the neighboring facet that shares `ridge` but does not include `apex`.
# ----------------------------------------------------------------------
def neighbor_across_ridge(X0, ridge, apex, interior):
    J = X0.shape[1]
    for d in range(J):
        if d == apex or d in ridge:
            continue
        S = sorted(ridge | {d})
        fn = facet_normal(X0, S, interior)
        if fn is None:
            continue
        if all_donors_in_halfspace(*fn, X0):
            return d, fn
    return None, None


def visible_facets(X0, X1):
    """
    All facets of conv(X0) visible from X1.
    """
    interior = interior_point(X0)
    seed = seed_visible_facet(X0, X1, interior)
    fn0 = facet_normal(X0, sorted(seed), interior)
    facets = {seed: fn0}
    q = deque([seed])
    while q:
        F = q.popleft()
        for apex in F:
            ridge = F - {apex}
            d, fn = neighbor_across_ridge(X0, ridge, apex, interior)
            if d is None:
                continue
            G = frozenset(ridge | {d})
            if G in facets:
                continue
            if visible(*fn, X1):
                facets[G] = fn
                q.append(G)
    return facets


# ----------------------------------------------------------------------
# Find the paper's allowed sign patterns on each face.
# Then keep the patterns whose closed residual regions are nonempty.
# ----------------------------------------------------------------------
LP_COUNT = [0]      # Number of LP solves, used in the face-count report.
RAY_TOL = 1e-9
_RAY_COMB_CAP = 2_000


def _normal_cone_yspace(face, X0):
    """
    Return (B, A, d) to describe the face's normal cone.

    B contains d independent directions perpendicular to the face.
    A @ y >= 0 ensures the normal c = B @ y puts every donor on or
    behind the plane through the face.
    """
    cols = sorted(face)
    K, J = X0.shape
    anchor = cols[0]
    edges = X0[:, cols[1:]] - X0[:, anchor, None]
    if edges.shape[1]:
        U, sv, _ = np.linalg.svd(edges)
        rank = int((sv > RAY_TOL * sv[0]).sum())
    else:
        U, rank = np.eye(K), 0
    B = U[:, rank:]
    others = [j for j in range(J) if j not in face]
    A = (X0[:, anchor, None] - X0[:, others]).T @ B
    return B, A, B.shape[1]


def _patterns_matching_normal(normal):
    """
    Each zero coordinate allows either sign.
    """
    choices = []
    for value in normal:
        if value > RAY_TOL:
            choices.append((1,))
        elif value < -RAY_TOL:
            choices.append((-1,))
        else:
            choices.append((-1, 1))
    return itertools.product(*choices)


def _normal_has_signs(B, A, signs):
    """
    Check for a nonzero normal c = B @ y with signs * c >= 0.
    """
    rows = np.vstack([A, signs[:, None] * B])
    # With matching signs, signs @ c is the sum of absolute coordinates.
    # Set it to one to exclude the zero normal.
    result = linprog(
        np.zeros(B.shape[1]),
        A_ub=-rows, b_ub=np.zeros(len(rows)),
        A_eq=(signs @ B)[None, :], b_eq=[1.0],
        bounds=(None, None), **LP,
    )
    LP_COUNT[0] += 1
    if result.status == 2:
        return False
    if result.status != 0:
        raise RuntimeError("Normal-cone LP failed: " + result.message)
    return True


def allowed_sign_patterns(face, X0, normal=None):
    """
    Return the paper's allowed patterns in {-1, +1}^K.

    Each pattern matches a nonzero normal to the face using weak signs.
    This step depends only on the donors, not on the treated unit.
    """
    K = X0.shape[0]
    if normal is not None and len(face) == K:
        length = np.linalg.norm(normal)
        if length == 0:
            return []
        return [np.array(s, dtype=int)
                for s in _patterns_matching_normal(normal / length)]

    B, A, d = _normal_cone_yspace(face, X0)
    rows = np.vstack([A, B])

    # For large enumerations, check the same condition by LP instead.
    if math.comb(len(rows), d - 1) > _RAY_COMB_CAP:
        return [
            np.array(s, dtype=int)
            for s in itertools.product((-1, 1), repeat=K)
            if _normal_has_signs(B, A, np.array(s))
        ]

    patterns = set()
    for selected in itertools.combinations(range(len(rows)), d - 1):
        equations = rows[list(selected)]
        _, sv, Vt = np.linalg.svd(equations)
        if sv.size and sv[-1] <= RAY_TOL * max(1.0, sv[0]):
            continue
        # The selected donor ties and coordinate zeros leave one direction.
        for y in (Vt[-1], -Vt[-1]):
            if np.any(A @ y < -RAY_TOL):
                continue
            patterns.update(_patterns_matching_normal(B @ y))
    return [np.array(s, dtype=int) for s in sorted(patterns)]


def piece_rows(signs, Xc, X1):
    """
    Return G, h for the weak residual constraints G @ weights <= h.
    """
    signs = np.asarray(signs)
    return signs[:, None] * Xc, signs * X1


def feasible_pieces(face, X0, X1, normal=None):
    """
    Keep allowed sign patterns with a nonempty closed region on the face.
    """
    Xc = X0[:, sorted(face)]
    m = Xc.shape[1]
    pieces = []
    for signs in allowed_sign_patterns(face, X0, normal):
        G, h = piece_rows(signs, Xc, X1)
        if np.any(G.min(axis=1) > h + TOL):
            continue
        result = linprog(
            np.zeros(m), A_ub=G, b_ub=h,
            A_eq=np.ones((1, m)), b_eq=[1.0], bounds=(0, None), **LP,
        )
        LP_COUNT[0] += 1
        if result.status == 2:
            continue
        if result.status != 0:
            raise RuntimeError("Residual-region LP failed: " + result.message)
        pieces.append(signs)
    return pieces


# ----------------------------------------------------------------------
# Branch and bound over visible faces.
# ----------------------------------------------------------------------
def _cannot_beat(bound, U):
    absolute_tolerance = 1e-9
    relative_tolerance = 1e-6
    tolerance = max(absolute_tolerance, relative_tolerance * abs(bound))
    return bound - tolerance >= U


def _finish(face, beta, J):
    w = np.zeros(J)
    w[sorted(face)] = beta
    w[w < 1e-7] = 0.0
    return w / w.sum()


def solve_exterior_hullfree(X0, X1, Z1, Z0, pieces_fn=feasible_pieces):
    K, J = X0.shape
    vfacets = visible_facets(X0, X1)
    facet_list = list(vfacets)

    def incident_visible(face):
        return [F for F in facet_list if face <= F]

    _, w_e = qp_simplex(X1, X0)
    w_e = np.clip(w_e, 0, None)
    w_e = w_e / w_e.sum()
    U = float((Z1 - Z0 @ w_e) @ (Z1 - Z0 @ w_e))
    best = (U, None, None)

    heap = []
    seen = set()
    tie = itertools.count()
    for F in facet_list:
        seen.add(F)
        Zc = Z0[:, sorted(F)]
        _, weights = qp_simplex(Z1, Zc)
        bound = simplex_lower_bound(Z1, Zc, weights)
        heapq.heappush(heap, (bound, next(tie), F))

    faces_touched = 0
    while heap:
        bound, _, face = heapq.heappop(heap)
        if _cannot_beat(bound, U):
            break
        faces_touched += 1
        fn = vfacets.get(face)                       # full facet => (n, b); else None
        normal = fn[0] if fn is not None else None
        pieces = pieces_fn(face, X0, X1, normal)
        cols = sorted(face)
        Zc, Xc = Z0[:, cols], X0[:, cols]
        for s in pieces:
            Gp, hp = piece_rows(s, Xc, X1)
            val, beta = qp_simplex(Z1, Zc, Gp, hp)
            if val is not None and val < U:
                U = val
                best = (val, face, beta)
        for v in face:
            sub = face - {v}
            if not sub or sub in seen or not incident_visible(sub):
                continue
            seen.add(sub)
            subcols = sorted(sub)
            if _cannot_beat(affine_bound(Z1, Z0[:, subcols]), U):
                continue
            Zc = Z0[:, subcols]
            _, weights = qp_simplex(Z1, Zc)
            bound = simplex_lower_bound(Z1, Zc, weights)
            if not _cannot_beat(bound, U):
                heapq.heappush(heap, (bound, next(tie), sub))

    loss, face, beta = best
    info = dict(
        faces_touched=faces_touched,
        n_visible_facets=len(vfacets),
        face=None if face is None else tuple(sorted(face)),
        visible_facets=[(tuple(sorted(F)), n, b) for F, (n, b) in vfacets.items()],
    )
    w_out = w_e if face is None else _finish(face, beta, J)
    # Check that removing numerical zero weights did not change the fit.
    realized = float((Z1 - Z0 @ w_out) @ (Z1 - Z0 @ w_out))
    if abs(realized - loss) > 1e-6 * max(1.0, abs(loss)):
        raise RuntimeError(
            "solve_exterior_hullfree internal inconsistency: reported loss "
            "%.10g but the returned weights realize %.10g. This indicates a "
            "solver bug; please report it."
            % (loss, realized))
    return realized, w_out, info


# ----------------------------------------------------------------------
# For general positive definite V, optimize over the closed visible facets.
# Use affine lower bounds to skip facets that cannot improve the current fit.
# ----------------------------------------------------------------------
def solve_exterior_pd(X0, X1, Z1, Z0):
    """
    Return the PD loss, donor weights, and search information.

    info["facet_vertices"] lists all donor indices in the selected facet,
    including those with zero weight. It is empty for the fallback fit.
    """
    K, J = X0.shape
    vfacets = visible_facets(X0, X1)

    heap = []
    tie = itertools.count()
    for F in vfacets:
        heapq.heappush(heap, (affine_bound(Z1, Z0[:, sorted(F)]), next(tie), F))

    U = float("inf")
    best = (U, None, None)
    facets_solved = 0
    while heap:
        bound, _, F = heapq.heappop(heap)
        if _cannot_beat(bound, U):
            break
        facets_solved += 1
        val, beta = qp_simplex(Z1, Z0[:, sorted(F)])
        if val is not None and val < U:
            U = val
            best = (val, F, beta)

    loss, F, beta = best
    info = dict(
        n_visible_facets=len(vfacets),
        facets_solved=facets_solved,
        V_class="pd",
        facet_vertices=sorted(F) if F is not None else [],
        face=None if F is None else tuple(sorted(F)),
        visible_facets=[(tuple(sorted(G)), n, b) for G, (n, b) in vfacets.items()],
    )
    if F is None:                                   # degenerate: Euclidean fallback
        _, w_e = qp_simplex(X1, X0)
        w_e = np.clip(w_e, 0, None)
        w_e = w_e / w_e.sum()
        return float((Z1 - Z0 @ w_e) @ (Z1 - Z0 @ w_e)), w_e, info
    weights = _finish(F, beta, J)
    residual = Z1 - Z0 @ weights
    return float(residual @ residual), weights, info
