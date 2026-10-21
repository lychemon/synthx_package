# `synthx`: exact synthetic control in Stata

Alberto Abadie and Yucheng Shang

Version: **0.9.1**.

Original code and documentation are licensed under the [MIT License](LICENSE).
Third-party datasets retain their original terms.

`synthx` implements the estimator in *Synthetic Control: An Exact Geometric
Formulation*. It minimizes validation loss over donor weights allowed by the
positive diagonal or full positive definite class of predictor weight matrices.
The exterior solver searches closed regions of visible faces, using the same
numerical routines as `synth-AS-replication`. Stata prepares the panel and displays
the results; Python performs the optimization.

## Requirements

Use Stata 16.1 or later and Python 3.12 or later. The tested dependencies are
NumPy 2.5.1, SciPy 1.18.0, CVXPY 1.9.3, and Clarabel 0.11.1. Figures use Stata,
so Matplotlib is not required.

Install the dependencies into the Python interpreter that Stata uses:

```bash
python -m pip install -r "<package>/code/python/requirements.txt"
```

Then configure Stata and check the installation:

```stata
python search
python set exec "C:/path/to/python.exe", permanently
adopath ++ "<package>/code/stata/bridge"
synthx_checkdeps
```

On macOS or Linux, use the path to your Python executable. `synthx_checkdeps`
reports package versions, imports the solver, and solves a small Clarabel QP.
It should finish with `PASS`. Use `python query` to inspect Stata's interpreter.

## Getting started

From `code/`, run the California example:

```stata
do example_california.do
```

It estimates the diagonal and PD models, reports predictor weights, draws the
outcome path, and runs a placebo test. For your own data, copy
`code/03_your_own_data.do` and replace the lines marked `<<<`.

The command takes a long panel with a numeric unit id and integer time variable:

```stata
xtset unit time
synthx outcome x1 x2 outcome(1988) outcome(1980), ///
    trunit(3) trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) ///
    vweights figure
```

| Option | Meaning |
|---|---|
| `trunit(#)` | Treated unit id. |
| `trperiod(#)` | First treated period; sets the default placebo post-period and figure line. It does not enter estimation. |
| `prewindow(lo hi)` | Validation outcomes used to choose donor weights. |
| `xperiod(lo hi)` | Window for averaging plain predictors; defaults to `prewindow()`. |
| `var(periods)` | Predictor averaged over those periods. For example, `beer(1984/1988)` or `schooling(1970 1975)`. |
| `unit()` and `time()` | Panel variables, if not supplied by `xtset`. |
| `vclass(diagonal)` | Optimize over positive diagonal metrics; the default. |
| `vclass(pd)` | Optimize over full positive definite metrics. |
| `vclass(invvar)` or `vclass(equal)` | Use a fixed inverse-variance or identity metric. |
| `vweights` | Report predictor-weight representatives and store their ranges. Exterior fits only. |
| `faces` | Count feasible, visible, and all faces by dimension. This checks every visible face and can be slow. Exterior fits only. |
| `allplacebos` | Refit each unit and rank the treated unit's post/pre-treatment MSPE ratio. |
| `postwindow(lo hi)` | Periods for the placebo numerator; defaults to `trperiod()` through the last panel period. |
| `figure` | Draw the treated and synthetic paths, or placebo gaps with `allplacebos`. |
| `treatline(#)` | Move the figure's vertical line; defaults to `trperiod()`. |
| `maxpre(#)` | Hide placebo paths with pre-treatment RMSE above this multiple of the treated RMSE. Does not change the test. |
| `title()`, `xtitle()`, `ytitle()`, `legend()` | Replace figure defaults. |
| `graphopts()` | Pass other Stata `twoway` options. |

For a leakage-free design, use predictor windows separate from validation:
Germany uses `xperiod(1971 1980) prewindow(1981 1990)`.

The placebo rank is one plus the number of ratios exceeding the treated ratio.
The p-value is that rank divided by the number of units that solved, including
the treated unit. Failed placebo fits are reported and excluded.

## Predictor weights and stored results

`vweights` describes the predictor matrices that reproduce the fitted donor
weights. Predictors are standardized across treated and donor units, and the
matrices have trace one. For each class it returns the log-determinant maximizer,
the analytic center, and the minimum and maximum of each diagonal entry.
A diagonal range of `[0, 0]` marks a predictor forced to zero in the closed set.
The center omits log barriers for donor constraints whose slack is always zero.

| Result | Contents |
|---|---|
| `e(b)` | Donor weights, with columns named by unit id. |
| `e(loss)`, `e(regime)`, `e(vclass)` | Validation loss, interior/exterior regime, and metric class. |
| `e(slack)` | Hull-test slack for a single fit. |
| `e(exact)` | Local geometry checks: 1 passed, 0 failed; absent when not applicable. This does not verify the full general-position assumption. |
| `e(balance)`, `e(plabels)` | Treated, synthetic, and donor-average predictors; labels in row order. |
| `e(Vdiag)`, `e(Vpd)` | K by 4: log-det representative, center, minimum, maximum. With `vweights`. |
| `e(pd_mineig)` | Smallest eigenvalue of the PD log-det representative. |
| `e(faces)` | K by 3: feasible / visible / all, rows for dimensions 0 through K-1. With `faces`. |
| `e(pvalue)`, `e(rank)`, `e(n)` | Placebo p-value, rank, and number of successful fits. |
| `e(placebos)` | Unit id, pre-RMSE, post-RMSE, and MSPE ratio. |
| `e(gaps)` | Period in the first column, followed by one gap path per successful unit. |

`help synthx` gives the full syntax. `code/synthx_tools.do` provides
`synthx_series` for outcome paths, `synthx_weights` for labeled donor weights,
and `synthx_unitid` for looking up a unit id. Built-in figures are named
`esc_path` and `esc_placebos`; use `graph export` to save them.

## Numerical methods and limits

The shared files are `simplex_solver.py`, `visible_facet_solver.py`, and
`predictor_weights.py`. QPs and predictor-weight problems use CVXPY/Clarabel;
LPs use SciPy/HiGHS. The diagonal exterior search uses binary signs and weak
residual inequalities, including boundary points in the closure. Lower bounds
allow branches to be skipped without solving every region.

The paper's general-position assumptions justify the exterior search.
`e(exact)` checks visible facets and the selected face only; passing these
checks is not a proof that the full assumption holds throughout the data.
Optimized classes need more donors than predictors. Missing predictors or
validation outcomes exclude a donor, with a count reported.

The diagonal class drops constant donor predictors with a note; PD retains them.
Nearly collinear predictors produce a note. The diagonal class stops on an
exactly lower-dimensional donor hull. The PD class reduces to the
donor span when the treated point is in that span; if it lies outside the span,
all donor-hull points are available and the solver minimizes loss over the full
simplex. A near-collinearity fallback to the inverse-variance metric is reported.
Fixed metrics retain their meaning on the original predictor scale.
Treated points exactly on the donor-hull boundary are not supported by the
current interior/exterior routing.

Small donor weights are set to zero and renormalized. The interior cell search
stops with an error if it exceeds 200,000 cells. Predictor-weight and face reports
are available for exterior fits; a report that cannot be computed is identified
in the Results window.

## Reproducing the paper

From `code/`, run `do run_all.do`. Outputs go to `output/`.

| Do-file | Outputs |
|---|---|
| `00_prepare_data.do` | `data/smoking.dta` and `data/germany.dta` from the CSV panels. |
| `01_california.do` | Tables 1-4 and D.1; Figures 5 and D.1; California results log and series. |
| `02_germany.do` | Tables 5-6 and D.2; Figures 6 and D.2; Germany results log and series. |
| `example_california.do` | Worked example, outcome-path figure, and placebo figure. |
| `03_your_own_data.do` | Template for another panel; runs on California until edited. |

Each do-file can also run separately from `code/`. The application scripts
check the paper values listed in their assertions. The original ADH weights
used for comparison are entered from the published papers.

The data are the California panel of Abadie, Diamond, and Hainmueller (2010),
39 states over 1970-2000, and their Germany panel (2015), 17 countries over
1960-2003. The preparation script creates alphabetical, value-labeled unit ids.

The Stata bridge is in `code/stata/bridge/`. The Python package is in
`code/python/synthx/` and can also be used directly with `from synthx import
estimate`; see `code/python/MANUAL.md`.
