# The exact synthetic control solver: Python API

Use Python 3.12 or later. From `code/python/`, install the pinned dependencies:

```bash
python -m pip install -r requirements.txt
```

## Estimate donor weights

```python
import sys
sys.path.insert(0, "<package>/code/python")
from synthx import estimate

result = estimate(X1, X0, Z1, Z0, V_class="diagonal")
weights = result["weights"]
loss = result["loss"]
```

`X1` has shape `(K,)`, `X0` has shape `(K, J)`, `Z1` has shape `(L,)`, and
`Z0` has shape `(L, J)`. Donor columns must have the same order in both matrices.

`estimate(X1, X0, Z1, Z0, V_class="diagonal", V=None)` supports:

- `V_class="diagonal"`: optimize over positive diagonal predictor metrics.
- `V_class="pd"`: optimize over full positive definite metrics.
- `V=`: use a fixed diagonal vector or positive definite matrix on the original
  predictor scale; the returned class is `"fixed"`.

The result contains `weights`, `loss`, `regime`, `V_class`, `interior_slack`,
and `info`. Weights are nonnegative and sum to one; loss is `||Z1 - Z0 @ w||²`.
`regime` is `"exterior"` or `"interior"`. `info` records search counts,
geometry diagnostics, and any fallback.

For exterior problems, the solver minimizes loss over closed regions of visible
faces using branch and bound. The diagonal regions have weak sign constraints;
the PD class uses the visible facets. QPs use CVXPY/Clarabel, and pruning uses
lower bounds from the shared `simplex_solver.py`. Interior problems use the
package's cell walks.

## Recover predictor weights

```python
from synthx.predictor_weights import diag_rationalizers, pd_rationalizers, count_faces

diagonal = diag_rationalizers(X1, X0, weights)
pd = pd_rationalizers(X1, X0, weights)
```

These functions standardize predictors across the treated and donor units and
normalize predictor matrices to trace one. They return full-precision values;
round only when displaying results. Set `standardized=True` when the inputs are
already standardized. `pd_rationalizers(..., ranges=False)` skips the range calculations.

| Function | Returned fields |
|---|---|
| `diag_rationalizers` | `logdet`, `center`, `lo`, `hi`, `support`, `implicit_eq`. |
| `pd_rationalizers` | `Vstar`, `center`, `logdet_diag`, `center_diag`, `lo`, `hi`, `min_eig`, `center_min_eig`, `implicit_eq`. |

`logdet` and `Vstar` maximize the log determinant on the feasible set.
`center` adds log barriers for inactive-donor slacks that can be positive.
Constraints whose slack is always zero are treated as equalities. `lo` and
`hi` give the range of each diagonal entry; a diagonal predictor with range
`[0, 0]` is forced to zero in the closed set. These calculations use the paper's
CVXPY formulations.

`count_faces(X1, X0, progress=None)` returns a K by 3 array. Rows represent
dimensions 0 through K-1; columns count faces with a feasible closed diagonal
region, faces below visible facets, and all hull faces. It assumes simplicial
faces and checks every visible face, independently of the optimization search.

## Assumptions and numerical limits

- The Stata interface requires more donors than predictors for optimized metrics.
- The exterior algorithm relies on the paper's general-position assumption.
  `info["general_position"]` checks visible facets and the selected face.
  Its `exact` flag reports those local checks, not a proof of the full assumption.
- The diagonal class drops constant donor predictors with a warning; PD retains
  them. Near collinearity produces a warning and can trigger a reported
  inverse-variance fallback.
- Exactly lower-dimensional donor hulls stop the diagonal solver. For PD,
  the solver reduces to the donor span if the treated point is in that span.
  If it is outside the span, the outcome QP uses all simplex weights.
- Treated points exactly on the donor-hull boundary are unsupported by the
  current interior/exterior routing.
- Small donor weights are set to zero and renormalized. Solver failures raise
  errors rather than returning an unchecked QP result.
- The interior search raises an error after 200,000 cells.

## Files

All solver modules are under `synthx/`; their relative imports keep them separate
from the paper's support modules when both are loaded in one Python process.

| File | Role |
|---|---|
| `sc_estimator.py` | Prepare inputs, choose the regime, and handle fixed metrics. |
| `visible_facet_solver.py` | Visible facets, closed sign regions, and exterior search. |
| `simplex_solver.py` | Clarabel simplex QPs and safe lower bounds. |
| `predictor_weights.py` | Predictor representatives, ranges, and face counts. |
| `interior_cellwalk.py`, `interior_vsearch.py` | Diagonal interior search and cell feasibility. |
| `interior_pd.py` | PD interior search. |
| `gpcheck.py` | Local general-position diagnostics. |
