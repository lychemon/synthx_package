"""
Check visible facets and the selected face for general-position violations.

These are local diagnostics, not a test of the full global assumption.
"""
from __future__ import annotations

import itertools

import numpy as np

REL_TOL = 1e-9


def _affinely_independent(P, rel_tol=REL_TOL):
    """
    P: t x m points. True if the t points are affinely independent.
    """
    t = P.shape[0]
    if t <= 1:
        return True
    M = np.column_stack([P, np.ones(t)])
    s = np.linalg.svd(M, compute_uv=False)
    return s.size >= t and s[-1] > rel_tol * max(s[0], 1.0)


def general_position_check(X0, X1, face, visible_facets, diagonal=True, labels=None, rel_tol=REL_TOL):
    """
    Check simplicial visible facets and independence on the selected face.

    For diagonal V, check each coordinate subset; for PD V, use all coordinates.
    Return exact, nonsimplicial, face_dependent, ties, and messages.
    """
    X0 = np.asarray(X0, float)
    K, J = X0.shape
    lab = (lambda j: str(labels[j])) if labels is not None else (lambda j: str(j))
    out = dict(exact=True, nonsimplicial=[], face_dependent=[], ties=[], messages=[])

    if not visible_facets or not face:
        out["exact"] = False
        out["messages"].append("Missing geometry for the local general-position checks.")

    visible_donors = set()
    for F, n, b in visible_facets or []:
        n = np.asarray(n, float).ravel()
        on = np.flatnonzero(np.abs(X0.T @ n - b) <= rel_tol * max(1.0, abs(float(b))) * max(1.0, float(np.abs(X0).max())))
        visible_donors.update(int(j) for j in F)
        if on.size > K:
            out["nonsimplicial"].append(sorted(int(j) for j in on))
    if out["nonsimplicial"]:
        out["exact"] = False
        out["messages"].append("%d visible facet(s) contain more than K donors (not simplicial): %s"
                               % (len(out["nonsimplicial"]), "; ".join("[" + ", ".join(lab(j) for j in f) + "]" for f in out["nonsimplicial"][:3])))

    if face:
        cols = [int(j) for j in face]
        P = X0[:, cols].T                         # donors of the face as rows
        subsets = ([tuple(range(K))] if not diagonal else
                   [S for m in range(1, K + 1) for S in itertools.combinations(range(K), m)])
        for S in subsets:
            m = len(S)
            tmax = min(m + 1, len(cols))
            for t in range(2, tmax + 1):
                for sub in itertools.combinations(range(len(cols)), t):
                    if not _affinely_independent(P[list(sub)][:, list(S)], rel_tol):
                        out["face_dependent"].append((list(S), [cols[i] for i in sub]))
            if len(out["face_dependent"]) > 20:
                break
        if out["face_dependent"]:
            out["exact"] = False
            S, d = out["face_dependent"][0]
            out["messages"].append("the donors of the solution face are not in general position: donors %s are "
                                   "affinely dependent on predictor(s) %s (%d such case(s))"
                                   % ([lab(j) for j in d], [k + 1 for k in S], len(out["face_dependent"])))

    vd = sorted(visible_donors) if visible_donors else list(range(J))
    for k in range(K):
        v = X0[k, vd]
        order = np.argsort(v)
        for a, b2 in zip(order[:-1], order[1:]):
            if abs(v[a] - v[b2]) <= rel_tol * max(1.0, abs(v[a])):
                out["ties"].append((k, vd[a], vd[b2]))
    return out
