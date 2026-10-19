# The exact synthetic control solver (Python)

This is the solver used by the Stata command `synthx`.

```python
from synthx import estimate

result = estimate(X1, X0, Z1, Z0, V_class="diagonal")  # or "pd"
weights = result["weights"]
loss = result["loss"]
```

`X1` and `X0` contain predictors, with one donor per column of `X0`.
`Z1` and `Z0` contain validation outcomes in the same donor order.
A fixed predictor metric can be supplied with `V=`.

The files `synthx/simplex_solver.py`, `synthx/visible_facet_solver.py`, and
`synthx/predictor_weights.py` follow the corresponding routines in
`synth-AS-replication/code/support/`. The package also includes the estimator,
interior solvers, and geometry diagnostics needed for general panel data.

Use Python 3.12 or later and install `requirements.txt`. The solver uses
NumPy, SciPy, CVXPY, and Clarabel. See [MANUAL.md](MANUAL.md) for the API and
numerical limits.

From the package root, run `python -m unittest discover -s tests` to check the
solver and bridge. When the sibling `synth-AS-replication` is present, this also
runs its solver regressions against the package implementation.
