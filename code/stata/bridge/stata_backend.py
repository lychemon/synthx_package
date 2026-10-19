"""
Prepare panel data, call the solver, and return results to Stata.

build() and fit() also work outside Stata. The other entry points use sfi to
store results as sx_* scalars, matrices, and locals for synthx.ado.
"""

from __future__ import annotations

import importlib
import os
import sys
import warnings

import numpy as np

# Load this installation's solver, including on repeated Stata calls.
_SOLVER_DIR = os.path.realpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
sys.path[:] = [p for p in sys.path if os.path.realpath(p) != _SOLVER_DIR]
sys.path.insert(0, _SOLVER_DIR)
for _name in list(sys.modules):
    if _name == "synthx" or _name.startswith("synthx."):
        del sys.modules[_name]
importlib.invalidate_caches()
for _name in ("synthx.sc_estimator", "synthx.predictor_weights"):
    importlib.import_module(_name)
for _name, _module in list(sys.modules.items()):
    if _name == "synthx" or _name.startswith("synthx."):
        _relative = _name.split(".")[1:]
        _expected = os.path.join(_SOLVER_DIR, "synthx", *_relative)
        _expected = (_expected + ".py") if _relative else os.path.join(_expected, "__init__.py")
        if os.path.normcase(os.path.realpath(_module.__file__)) != os.path.normcase(os.path.realpath(_expected)):
            raise ImportError("synthx: loaded %s from an unexpected path: %s" % (_name, _module.__file__))

from synthx.sc_estimator import estimate   # noqa: E402

_LASTFIT = None   # X1 / X0 / w / regime of the last treated-unit fit (for extras())
_PCACHE = None    # panel and per-unit outcome series cached for the placebo loop


# Plain Python: data preparation and estimation
def parse_special(s):
    """
    'var y1 y2 ...; var y; ...' -> [(var, [years]), ...].
    """
    out = []
    for part in (s or "").split(";"):
        toks = part.split()
        if len(toks) < 2:
            continue
        out.append((toks[0], [int(t) for t in toks[1:]]))
    return out


def predictor_labels(reg, specials):
    """
    Predictor labels in the row order build() produces: plain predictors,
    then specials with their periods, e.g. cigsale(1988) or beer(1984-1988).
    """
    labels = list(reg)
    for (v, yrs) in specials:
        if len(yrs) == 1:
            labels.append("%s(%d)" % (v, yrs[0]))
        elif yrs == list(range(yrs[0], yrs[-1] + 1)):
            labels.append("%s(%d-%d)" % (v, yrs[0], yrs[-1]))
        else:
            labels.append("%s(%s)" % (v, ",".join(str(y) for y in yrs)))
    return labels


def build(data, depvar, reg, specials, unitvar, timevar, trunit, plo, phi, xlo, xhi):
    """
    Predictor and validation matrices from a long panel.

    data: dict varname -> length-N float array (missing as NaN). Returns
    dict(donors, X1, X0, Z1, Z0, n_dropped): donors are the unit ids with
    complete predictors and validation outcomes, in increasing order; n_dropped
    counts the units (other than the treated one) dropped for missing values.
    """
    unit = data[unitvar]; time = data[timevar]; y = data[depvar]
    units = np.unique(unit)
    if not np.any(units == trunit):
        raise ValueError("trunit(%g) is not a value of the unit variable %s"
                         % (trunit, unitvar))
    years = np.arange(int(plo), int(phi) + 1)
    K = len(reg) + len(specials)
    Xall = np.full((K, len(units)), np.nan)
    Zall = np.full((len(years), len(units)), np.nan)

    def wmean(col, mask):
        v = col[mask & ~np.isnan(col)]
        return v.mean() if v.size else np.nan

    for i, u in enumerate(units):
        mu = unit == u
        r = 0
        for v in reg:
            Xall[r, i] = wmean(data[v], mu & (time >= xlo) & (time <= xhi)); r += 1
        for (v, yrs) in specials:
            Xall[r, i] = wmean(data[v], mu & np.isin(time, yrs)); r += 1
        for j, yr in enumerate(years):
            Zall[j, i] = wmean(y, mu & (time == yr))

    tc = int(np.where(units == trunit)[0][0])
    valid = ~(np.isnan(Xall).any(axis=0) | np.isnan(Zall).any(axis=0))
    if not valid[tc]:
        raise ValueError("the treated unit has missing predictor or outcome values "
                         "in its windows")
    don = np.array([i for i in range(len(units)) if i != tc and valid[i]])
    return dict(donors=units[don], X1=Xall[:, tc], X0=Xall[:, don],
                Z1=Zall[:, tc], Z0=Zall[:, don],
                n_dropped=int((~valid).sum()))


def _estimate_vclass(d, vclass):
    """
    Optimize over diagonal or PD metrics, or use a fixed metric.
    """
    K, J = np.asarray(d["X0"]).shape
    if vclass in ("diagonal", "pd") and J <= K:
        raise ValueError("the optimized-V classes need more donors than predictors "
                         "(J = %d donors with complete data, K = %d predictors); use "
                         "vclass(invvar) or vclass(equal), or drop predictors" % (J, K))
    if vclass == "equal":
        return estimate(d["X1"], d["X0"], d["Z1"], d["Z0"], V=np.ones(K))
    if vclass == "invvar":
        pooled = np.column_stack([d["X0"], np.asarray(d["X1"]).reshape(-1, 1)])
        v = np.var(pooled, axis=1)
        v[v == 0] = 1.0
        return estimate(d["X1"], d["X0"], d["Z1"], d["Z0"], V=1.0 / v)
    return estimate(d["X1"], d["X0"], d["Z1"], d["Z0"], V_class=vclass,
                    donor_labels=[("%d" % u) if float(u).is_integer() else ("%g" % u) for u in d["donors"]])


def fit(data, depvar, reg, specials, unitvar, timevar, trunit, plo, phi, xlo, xhi,
        vclass="diagonal"):
    """
    Fit the treated unit and return results with any solver warnings.
    """
    d = build(data, depvar, reg, specials, unitvar, timevar, trunit, plo, phi, xlo, xhi)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = _estimate_vclass(d, vclass)
    return dict(donors=d["donors"], weights=r["weights"], loss=r["loss"],
                regime=r["regime"], slack=r["interior_slack"],
                n_dropped=d["n_dropped"], X1=d["X1"], X0=d["X0"],
                notes=[str(c.message) for c in caught],
                gp=r["info"].get("general_position"))


def unit_series(data, unitvar, timevar, depvar, years):
    """
    dict unit -> outcome path over `years` (NaN where missing).
    """
    unit = data[unitvar]; time = data[timevar]; y = data[depvar]
    out = {}
    for u in np.unique(unit):
        s = np.full(len(years), np.nan)
        mu = unit == u
        for j, yr in enumerate(years):
            v = y[mu & (time == yr)]
            v = v[~np.isnan(v)]
            s[j] = v.mean() if v.size else np.nan
        out[u] = s
    return out


def placebo_row(data, depvar, reg, specials, unitvar, timevar, u,
                plo, phi, xlo, xhi, years, Y, post_lo, post_hi, vclass):
    """
    Fit unit u as if treated; return its gap path and pre/post RMSE.
    """
    d = build(data, depvar, reg, specials, unitvar, timevar, u, plo, phi, xlo, xhi)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = _estimate_vclass(d, vclass)
    synth = np.zeros(len(years))
    for k, dn in enumerate(d["donors"]):
        synth = synth + r["weights"][k] * Y[dn]
    gap = Y[u] - synth
    pre = np.sqrt(np.nanmean(gap[(years >= plo) & (years <= phi)] ** 2))
    post = np.sqrt(np.nanmean(gap[(years >= post_lo) & (years <= post_hi)] ** 2))
    ratio = (post ** 2) / (pre ** 2) if pre > 0 else np.inf
    return dict(unit=u, regime=r["regime"], pre=pre, post=post, ratio=ratio,
                gap=gap, weights=r["weights"], donors=d["donors"],
                loss=r["loss"], X1=d["X1"], X0=d["X0"],
                notes=[str(c.message) for c in caught],
                gp=r["info"].get("general_position"))


# Stata side (sfi). Arguments arrive as strings.
def _say(line):
    from sfi import SFIToolkit
    SFIToolkit.displayln(line)


def _notes(notes, prefix="synthx: note: "):
    for n in notes:
        _say(prefix + " ".join(str(n).split()))


def _fail(msg):
    """
    Report an error to Stata without a Python traceback: the .ado checks
    sx_err after each call and exits with the message.
    """
    from sfi import Scalar, Macro
    Scalar.setValue("sx_err", 1.0)
    Macro.setLocal("sx_errmsg", " ".join(str(msg).split())[:240])


def _read_panel(depvar, reg, specials, unitvar, timevar):
    from sfi import Data
    needed = list(dict.fromkeys([unitvar, timevar, depvar] + reg + [s[0] for s in specials]))
    M = np.array(Data.get(needed), dtype=float)
    M[np.abs(M) > 1e300] = np.nan          # Stata missing -> NaN
    return {nm: M[:, k] for k, nm in enumerate(needed)}


def _store_fit(res, reg, specials):
    """
    Store the treated fit and local geometry checks in Stata.
    """
    from sfi import Matrix, Scalar, Macro
    global _LASTFIT
    w = np.asarray(res["weights"], float).reshape(1, -1)
    Matrix.create("sx_W", 1, w.shape[1], 0)
    Matrix.store("sx_W", w)
    Matrix.setColNames("sx_W", [str(int(u)) for u in res["donors"]])
    Scalar.setValue("sx_loss", float(res["loss"]))
    Macro.setLocal("sx_regime", res["regime"])
    X1 = np.asarray(res["X1"], float); X0 = np.asarray(res["X0"], float)
    BAL = np.column_stack([X1, X0 @ w.ravel(), X0.mean(axis=1)])
    Matrix.create("sx_BAL", BAL.shape[0], 3, 0)
    Matrix.store("sx_BAL", BAL)
    Matrix.setColNames("sx_BAL", ["treated", "synthetic", "donoravg"])
    Macro.setLocal("sx_plabels", ";".join(predictor_labels(reg, specials)))
    gp = res.get("gp")
    if gp is None:                      # interior or fixed-V fit: not applicable
        Scalar.setValue("sx_exact", -1.0)
        Macro.setLocal("sx_gpnote", "")
    else:
        Scalar.setValue("sx_exact", 1.0 if gp["exact"] else 0.0)
        Macro.setLocal("sx_gpnote", " | ".join(gp["messages"])[:2000])
    _LASTFIT = dict(X1=X1, X0=X0, w=w.ravel(), regime=res["regime"])


def run(depvar, regpreds, special, unitvar, timevar,
        trunit, trperiod, plo, phi, xlo, xhi, vclass="diagonal", post_lo="", post_hi=""):
    """
    Single estimate. Stores sx_W, sx_loss, sx_slack, sx_regime, sx_BAL,
    sx_plabels, sx_ndropped and the path matrix sx_PATH (year, treated,
    synthetic, gap) for the figure. post_lo/post_hi are accepted for a uniform
    signature with the placebo calls and are not used here.
    """
    from sfi import Matrix, Scalar, Macro
    Scalar.setValue("sx_err", 0.0)
    try:
        reg = regpreds.split()
        specials = parse_special(special)
        trunit = float(trunit); plo = int(plo); phi = int(phi)
        xlo = int(xlo) if xlo else plo
        xhi = int(xhi) if xhi else phi
        data = _read_panel(depvar, reg, specials, unitvar, timevar)
        res = fit(data, depvar, reg, specials, unitvar, timevar,
                  trunit, plo, phi, xlo, xhi, vclass=(vclass or "diagonal"))
    except Exception as e:
        _fail(e)
        return
    _notes(res["notes"])
    _store_fit(res, reg, specials)
    Scalar.setValue("sx_slack", float(res["slack"]))
    Macro.setLocal("sx_ndropped", str(res["n_dropped"]))

    years = np.unique(data[timevar][data[timevar] >= plo])
    Y = unit_series(data, unitvar, timevar, depvar, years)
    synth = np.zeros(len(years))
    for k, dn in enumerate(res["donors"]):
        synth = synth + res["weights"][k] * Y[dn]
    PATH = np.column_stack([years, Y[trunit], synth, Y[trunit] - synth])
    Matrix.create("sx_PATH", PATH.shape[0], 4, 0)
    Matrix.store("sx_PATH", PATH)
    Matrix.setColNames("sx_PATH", ["year", "treated", "synth", "gap"])


def placebo_init(depvar, regpreds, special, unitvar, timevar,
                 trunit, trperiod, plo, phi, xlo, xhi, vclass="diagonal",
                 post_lo="", post_hi=""):
    """
    Read the panel and build the per-unit outcome series once; the .ado then
    calls placebo_one() for every unit so the Results window updates live.
    """
    global _PCACHE
    from sfi import Scalar
    Scalar.setValue("sx_err", 0.0)
    try:
        reg = regpreds.split()
        specials = parse_special(special)
        trunit_f = float(trunit); trperiod_i = int(trperiod)
        plo = int(plo); phi = int(phi)
        xlo = int(xlo) if xlo else plo
        xhi = int(xhi) if xhi else phi
        data = _read_panel(depvar, reg, specials, unitvar, timevar)
        if not np.any(np.unique(data[unitvar]) == trunit_f):
            raise ValueError("trunit(%g) is not a value of the unit variable %s" % (trunit_f, unitvar))
        years = np.unique(data[timevar][data[timevar] >= plo])
        _PCACHE = dict(data=data, depvar=depvar, reg=reg, specials=specials,
                       unitvar=unitvar, timevar=timevar, trunit=trunit_f,
                       plo=plo, phi=phi, xlo=xlo, xhi=xhi, years=years,
                       Y=unit_series(data, unitvar, timevar, depvar, years),
                       post_lo=int(post_lo) if post_lo else trperiod_i,
                       post_hi=int(post_hi) if post_hi else int(years.max()),
                       vclass=(vclass or "diagonal"))
    except Exception as e:
        _fail(e)


def placebo_one(oneunit, oneidx="", onetot=""):
    """
    Fit one placebo unit from the cache. Stores sx_one_ok (1/0), sx_one_pre,
    sx_one_post, sx_one_ratio, sx_one_gap, sx_one_years and, for the treated
    unit, the same results as run(). Prints a progress bar: one mark per unit.
    """
    from sfi import Matrix, Scalar, SFIToolkit
    if onetot:
        try:
            ntot = int(float(onetot)); ki = int(float(oneidx))
        except Exception:
            ntot = ki = 0
        if str(oneidx) == "1":
            SFIToolkit.displayln("")
            SFIToolkit.display("synthx: %d refits (%d placebos)  " % (ntot, max(ntot - 1, 0)))
        SFIToolkit.display(str(ki) if (ntot and (ki % 5 == 0 or ki == ntot)) else ".")
    if _PCACHE is None:
        Scalar.setValue("sx_one_ok", 0.0)
        return
    C = _PCACHE
    u = float(oneunit)
    try:
        row = placebo_row(C["data"], C["depvar"], C["reg"], C["specials"],
                          C["unitvar"], C["timevar"], u, C["plo"], C["phi"],
                          C["xlo"], C["xhi"], C["years"], C["Y"],
                          C["post_lo"], C["post_hi"], C["vclass"])
    except Exception as e:
        Scalar.setValue("sx_one_ok", 0.0)
        Scalar.setValue("sx_one_err", 1.0)
        from sfi import Macro
        Macro.setLocal("sx_one_errmsg", " ".join(str(e).split())[:200])
        return
    Scalar.setValue("sx_one_ok", 1.0)
    Scalar.setValue("sx_one_err", 0.0)
    Scalar.setValue("sx_one_pre", float(row["pre"]))
    Scalar.setValue("sx_one_post", float(row["post"]))
    rr = float(row["ratio"])
    Scalar.setValue("sx_one_ratio", rr if np.isfinite(rr) else 1e300)
    gap = np.asarray(row["gap"], float).reshape(1, -1)
    Matrix.create("sx_one_gap", 1, gap.shape[1], 0); Matrix.store("sx_one_gap", gap)
    yr = np.asarray(C["years"], float).reshape(1, -1)
    Matrix.create("sx_one_years", 1, yr.shape[1], 0); Matrix.store("sx_one_years", yr)
    if row["notes"]:
        SFIToolkit.displayln("")
        _notes(row["notes"], prefix="synthx: note (unit %g): " % u)
    if int(round(u)) == int(round(C["trunit"])):
        _store_fit(row, C["reg"], C["specials"])


# Predictor weights and face counts for the last treated fit.
def extras(vweights="", faces=""):
    """
    vweights -> sx_VDIAG (K x 4: log-det maximizer, analytic center, min, max
    over the diagonal rationalizers) and sx_VPD (K x 4: the same four columns
    over the positive definite rationalizers, diagonals trace-normalized) plus
    sx_pd_mineig; faces -> sx_FACES (K x 3, rows = dimension 0..K-1: attainable
    under diagonal V / visible / all). Flags sx_vdiag_ok, sx_vpd_ok, sx_faces_ok.
    """
    from sfi import Matrix, Scalar
    for nm in ("sx_vdiag_ok", "sx_vpd_ok", "sx_faces_ok"):
        Scalar.setValue(nm, 0.0)
    if _LASTFIT is None:
        return
    X1, X0, w = _LASTFIT["X1"], _LASTFIT["X0"], _LASTFIT["w"]
    if _LASTFIT["regime"] != "exterior":
        _say("synthx: vweights and faces apply to the exterior regime only "
             "(the treated unit is interior); skipped.")
        return
    if vweights:
        from synthx.predictor_weights import diag_rationalizers, pd_rationalizers
        try:
            dg = diag_rationalizers(X1, X0, w)
            VD = np.column_stack([dg["logdet"], dg["center"], dg["lo"], dg["hi"]])
            Matrix.create("sx_VDIAG", VD.shape[0], 4, 0)
            Matrix.store("sx_VDIAG", VD)
            Matrix.setColNames("sx_VDIAG", ["logdet", "center", "min", "max"])
            Scalar.setValue("sx_vdiag_ok", 1.0)
        except Exception as e:
            _say("synthx: diagonal predictor weights not computed (%s)" % e)
        try:
            pd_ = pd_rationalizers(X1, X0, w)
            VP = np.column_stack([pd_["logdet_diag"], pd_["center_diag"], pd_["lo"], pd_["hi"]])
            Matrix.create("sx_VPD", VP.shape[0], 4, 0)
            Matrix.store("sx_VPD", VP)
            Matrix.setColNames("sx_VPD", ["logdet", "center", "min", "max"])
            Scalar.setValue("sx_pd_mineig", float(pd_["min_eig"]))
            Scalar.setValue("sx_vpd_ok", 1.0)
        except Exception as e:
            _say("synthx: positive definite predictor weights not computed (%s)" % e)
    if faces:
        from synthx.predictor_weights import count_faces
        _say("synthx: counting the faces of the donor hull (attainable / visible / all) ...")

        def _prog(i, n):
            _say("synthx:   %d / %d visible faces tested" % (i, n))
        try:
            F = count_faces(X1, X0, progress=_prog)
            Matrix.create("sx_FACES", F.shape[0], 3, 0)
            Matrix.store("sx_FACES", np.asarray(F, float))
            Matrix.setColNames("sx_FACES", ["attainable", "visible", "all"])
            Scalar.setValue("sx_faces_ok", 1.0)
        except Exception as e:
            _say("synthx: face counts not computed (%s)" % e)
