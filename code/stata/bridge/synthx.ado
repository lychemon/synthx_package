*! version 0.9.0  19oct2026
*! synthx -- exact synthetic control (Abadie and Shang, Synthetic Control: An Exact Geometric Formulation).
*! Stata prepares the data and calls the Python solver through stata_backend.py.
*! Run synthx_checkdeps to check Python and its dependencies; see help synthx for syntax.

* Allow the application do-files to load the command again.
capture program drop synthx
program synthx, eclass
    version 16.1
    // predictors (with inline special var(periods)) come before the comma;
    // syntax varlist would reject the parentheses, so parse manually as synth does.
    gettoken predspec optspec : 0, parse(",")
    local 0 `optspec'
    syntax , TRUnit(integer) TRPeriod(integer) ///
        PREwindow(numlist min=2 max=2 integer ascending) ///
        [ XPeriod(numlist min=2 max=2 integer ascending) ///
          POSTwindow(numlist min=2 max=2 integer ascending) ///
          Unit(varname) Time(varname) ALLPlacebos Figure Vclass(string) ///
          VWeights FACes TREATline(numlist max=1 integer) MAXPre(numlist max=1 >0) ///
          title(string asis) xtitle(string asis) ytitle(string asis) legend(string asis) ///
          graphopts(string asis) ]

    if "`vclass'" == "" local vclass "diagonal"
    if "`treatline'" == "" local treatline `trperiod'   // vertical line of the figure
    if !inlist("`vclass'", "diagonal", "pd", "invvar", "equal") {
        di as error "vclass() must be diagonal, pd, invvar, or equal"
        exit 198
    }

    gettoken depvar preds : predspec
    confirm numeric variable `depvar'
    local regpreds ""
    local special ""
    while `"`preds'"' != "" {
        gettoken tok preds : preds, bind
        if strpos(`"`tok'"', "(") {
            local vn = substr(`"`tok'"', 1, strpos(`"`tok'"', "(") - 1)
            local in = substr(`"`tok'"', strpos(`"`tok'"', "(") + 1, .)
            local in = subinstr(`"`in'"', ")", "", .)
            numlist "`in'"
            local special `"`special'`vn' `r(numlist)'; "'
        }
        else {
            local regpreds `regpreds' `tok'
        }
    }
    if "`regpreds'" == "" & `"`special'"' == "" {
        di as error "specify at least one predictor"
        exit 198
    }
    if "`unit'" == "" | "`time'" == "" {
        qui xtset
        if "`unit'" == "" local unit = r(panelvar)
        if "`time'" == "" local time = r(timevar)
    }
    tokenize "`prewindow'"
    local plo `1'
    local phi `2'
    local xlo : word 1 of `xperiod'
    local xhi : word 2 of `xperiod'
    // The placebo post-period defaults to trperiod through the last panel year.
    local pslo : word 1 of `postwindow'
    local pshi : word 2 of `postwindow'
    // Locate the bridge using the global path or adopath.
    if "$synthx_bridge" != "" {
        local pp "$synthx_bridge"
    }
    else {
        capture qui findfile stata_backend.py
        if _rc {
            di as error "synthx: cannot find the bridge folder (stata_backend.py). Either"
            di as error `"  adopath ++ "<package>/code/stata/bridge"  or  global synthx_bridge "<that folder>""'
            exit 601
        }
        local pp = subinstr("`r(fn)'", "stata_backend.py", "", .)
    }
    // a trailing separator would end the raw Python string with a backslash
    while inlist(substr("`pp'", -1, 1), char(92), "/") {
        local pp = substr("`pp'", 1, length("`pp'") - 1)
    }

    // Check Python and package presence. synthx_checkdeps also tests versions and Clarabel.
    tempname pfok
    capture python: import sys, importlib.util as _u; from sfi import Scalar; Scalar.setValue("`pfok'", int(sys.version_info >= (3, 12) and all(_u.find_spec(p) is not None for p in ("numpy", "scipy", "cvxpy", "clarabel"))))
    if _rc {
        di as error "synthx: no usable Python is configured in Stata. " ///
            "Run {bf:synthx_checkdeps} for setup help."
        exit 198
    }
    if scalar(`pfok') == 0 {
        di as error "synthx: Python 3.12+ with numpy, scipy, cvxpy, and clarabel is required. " ///
            "Run {bf:synthx_checkdeps} for details."
        exit 198
    }

    // Load this bridge file afresh; later calls share its state within this command.
    python: import os, importlib.util; _spec = importlib.util.spec_from_file_location("_synthx_stata_backend", os.path.join(r"`pp'", "stata_backend.py")); _sb = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_sb)

    if "`allplacebos'" == "" {
        // Fit the treated unit.
        python: from sfi import Macro; _sb.run(Macro.getLocal("depvar"), Macro.getLocal("regpreds"), Macro.getLocal("special"), Macro.getLocal("unit"), Macro.getLocal("time"), Macro.getLocal("trunit"), Macro.getLocal("trperiod"), Macro.getLocal("plo"), Macro.getLocal("phi"), Macro.getLocal("xlo"), Macro.getLocal("xhi"), Macro.getLocal("vclass"), Macro.getLocal("pslo"), Macro.getLocal("pshi"))
        if scalar(sx_err) == 1 {
            di as error "synthx: `sx_errmsg'"
            capture scalar drop sx_err
            exit 498
        }
    }
    else {
        // Fit one unit per Python call so Stata can display progress.
        quietly levelsof `unit', local(_ulist)
        local _n : word count `_ulist'
        local _k = 0
        local _m = 0
        local _ny = 0
        local _trrow = 0
        local _nskip = 0
        // Cache the panel and outcome series for all placebo fits.
        python: from sfi import Macro; _sb.placebo_init(Macro.getLocal("depvar"), Macro.getLocal("regpreds"), Macro.getLocal("special"), Macro.getLocal("unit"), Macro.getLocal("time"), Macro.getLocal("trunit"), Macro.getLocal("trperiod"), Macro.getLocal("plo"), Macro.getLocal("phi"), Macro.getLocal("xlo"), Macro.getLocal("xhi"), Macro.getLocal("vclass"), Macro.getLocal("pslo"), Macro.getLocal("pshi"))
        if scalar(sx_err) == 1 {
            di as error "synthx: `sx_errmsg'"
            capture scalar drop sx_err
            exit 498
        }
        // Print progress inside Python to avoid duplicate lines in Stata.
        foreach _u of local _ulist {
            local ++_k
            local oneunit `_u'
            local oneidx  `_k'
            local onetot  `_n'
            python: from sfi import Macro; _sb.placebo_one(Macro.getLocal("oneunit"), Macro.getLocal("oneidx"), Macro.getLocal("onetot"))
            if scalar(sx_one_ok) == 0 {
                di as txt _n "synthx:     unit `_u' skipped: `sx_one_errmsg'"
                local ++_nskip
                continue
            }
            local ++_m
            if `_m' == 1 {
                local _ny = colsof(sx_one_years)
                matrix _INF = J(`_n', 4, .)
                matrix _GAP = J(`_ny', `=`_n'+1', .)
                forvalues t = 1/`_ny' {
                    matrix _GAP[`t', 1] = sx_one_years[1, `t']
                }
            }
            matrix _INF[`_m', 1] = `_u'
            matrix _INF[`_m', 2] = scalar(sx_one_pre)
            matrix _INF[`_m', 3] = scalar(sx_one_post)
            matrix _INF[`_m', 4] = scalar(sx_one_ratio)
            forvalues t = 1/`_ny' {
                matrix _GAP[`t', `=`_m'+1'] = sx_one_gap[1, `t']
            }
            if `_u' == `trunit' {
                local _trrow = `_m'
                matrix _Wtr   = sx_W
                matrix _BALtr = sx_BAL
                local _losstr = scalar(sx_loss)
                local _regtr  "`sx_regime'"
                local _platr  "`sx_plabels'"
            }
        }
        if `_trrow' == 0 {
            di as error "synthx: the treated unit could not be fit in the permutation"
            exit 498
        }
        // Keep the units that solved and assemble the stored results.
        if `_m' < `_n' {
            matrix _INF = _INF[1..`_m', 1..4]
            matrix _GAP = _GAP[1..`_ny', 1..`=`_m'+1']
        }
        if `_nskip' > 0 di as txt _n "synthx: `_nskip' unit(s) skipped; the p-value denominator is the `_m' units that solved."
        local _greater = 0
        forvalues i = 1/`_m' {
            if _INF[`i', 4] > _INF[`_trrow', 4] local ++_greater
        }
        matrix  sx_INF   = _INF
        matrix  sx_GAP   = _GAP
        matrix  sx_W     = _Wtr
        matrix  sx_BAL   = _BALtr
        scalar  sx_loss  = `_losstr'
        scalar  sx_rank  = `_greater' + 1
        scalar  sx_n     = `_m'
        scalar  sx_pval  = (`_greater' + 1) / `_m'
        scalar  sx_trcol = `_trrow'
        local   sx_regime  "`_regtr'"
        local   sx_plabels "`_platr'"
        capture matrix drop _INF _GAP _Wtr _BALtr
        capture scalar drop sx_one_ok sx_one_pre sx_one_post sx_one_ratio
        capture matrix drop sx_one_gap sx_one_years
    }

    // results placed by the backend.  Capture before ereturn post clears r().
    tempname b bal
    matrix `b' = sx_W
    matrix `bal' = sx_BAL
    local plabels "`sx_plabels'"
    local loss   = scalar(sx_loss)
    local regime "`sx_regime'"
    capture local slack = scalar(sx_slack)

    // Add predictor weights and face counts when requested.
    if "`vweights'`faces'" != "" {
        python: _sb.extras("`vweights'", "`faces'")
    }

    di as txt _n "Exact synthetic control (" as res "`regime'" as txt " case)"
    di as txt "Treated unit " as res `trunit' as txt ", validation loss " ///
       as res %9.4f `loss'

    // Display predictor balance.
    di as txt _n "Predictor balance:"
    di as txt _col(39) "Treated" _col(51) "Synthetic" _col(64) "Donor avg."
    di as txt "{hline 73}"
    local K = rowsof(`bal')
    local rest "`plabels'"
    forvalues i = 1/`K' {
        local pos = strpos("`rest'", ";")
        if `pos' {
            local lab = substr("`rest'", 1, `pos' - 1)
            local rest = substr("`rest'", `pos' + 1, .)
        }
        else local lab "`rest'"
        local t  = `bal'[`i', 1]
        local s  = `bal'[`i', 2]
        local da = `bal'[`i', 3]
        di as txt _col(5) "`lab'" as res ///
           _col(34) %12.3f `t' _col(48) %12.3f `s' _col(62) %12.3f `da'
    }
    di as txt "{hline 73}"
    if "`sx_ndropped'" != "" & "`sx_ndropped'" != "0" {
        di as txt "(`sx_ndropped' unit(s) dropped from the donor pool: missing predictor or outcome values)"
    }
    // Display local geometry diagnostics.
    capture local exact = scalar(sx_exact)
    if "`exact'" == "" local exact -1
    if `exact' == 0 {
        di as txt _n "synthx: note: the local geometry checks did not pass (e(exact) = 0)."
        di as txt "  These diagnostics do not verify the paper's full general-position assumption."
        di as txt "  `sx_gpnote'"
    }

    // Display nonzero donor weights, largest first.
    di as txt _n "Donor weights (nonzero):"
    local names : colnames `b'
    local J = colsof(`b')
    local ulab : value label `unit'
    local left ""
    forvalues j = 1/`J' {
        if abs(`b'[1, `j']) > 1e-6 local left `left' `j'
    }
    while "`left'" != "" {
        local jmax ""
        local wmax -1
        foreach j of local left {
            if `b'[1, `j'] > `wmax' {
                local wmax = `b'[1, `j']
                local jmax `j'
            }
        }
        local left : list left - jmax
        local nm : word `jmax' of `names'
        local lab ""
        if "`ulab'" != "" local lab : label `ulab' `nm'
        if "`lab'" == "" | "`lab'" == "`nm'" {
            di as txt "  " as res %-28s "unit `nm'" _col(44) as res %7.4f `wmax'
        }
        else {
            di as txt "  " as res %-28s "`lab'" as txt " (id " as res "`nm'" as txt ")" _col(44) as res %7.4f `wmax'
        }
    }

    // Display predictor weights.
    local showv 0
    if "`vweights'" != "" {
        if scalar(sx_vdiag_ok) == 1 | scalar(sx_vpd_ok) == 1 local showv 1
    }
    if `showv' {
        di as txt _n "Predictor weights that rationalize the donor weights"
        di as txt "(standardized predictors, trace one; each column sums to one)"
        di as txt _col(31) "Diagonal class" _col(53) "Positive definite class"
        di as txt _col(31) "log-det   center" _col(53) "log-det   center"
        di as txt "{hline 70}"
        local rest "`plabels'"
        forvalues i = 1/`K' {
            local pos = strpos("`rest'", ";")
            if `pos' {
                local lab = substr("`rest'", 1, `pos' - 1)
                local rest = substr("`rest'", `pos' + 1, .)
            }
            else local lab "`rest'"
            di as txt _col(3) substr("`lab'", 1, 26) _continue
            if scalar(sx_vdiag_ok) == 1 {
                di as res _col(31) %7.3f sx_VDIAG[`i', 1] _col(41) %7.3f sx_VDIAG[`i', 2] _continue
            }
            else di as txt _col(31) "    n/a" _continue
            if scalar(sx_vpd_ok) == 1 {
                di as res _col(53) %7.3f sx_VPD[`i', 1] _col(63) %7.3f sx_VPD[`i', 2]
            }
            else di as txt _col(53) "    n/a"
        }
        di as txt "{hline 70}"
        if scalar(sx_vdiag_ok) == 1 {
            di as txt "Diagonal: log-det maximizer and analytic center over the predictors not forced"
            di as txt "  to zero (Section 4.5.2 of the paper); a forced zero shows as 0 in both columns."
        }
        if scalar(sx_vpd_ok) == 1 {
            di as txt "Positive definite: trace-normalized diagonals of the log-det maximizer (smallest"
            di as txt "  eigenvalue " as res %8.2e scalar(sx_pd_mineig) as txt ") and of the analytic center."
        }
        di as txt "The identified range of each entry is stored in columns 3-4 of e(Vdiag) and e(Vpd)."
    }

    // Display face counts by dimension.
    if "`faces'" != "" {
        if scalar(sx_faces_ok) == 1 {
            di as txt _n "Faces of the donor hull, by dimension:"
            di as txt _col(5) "dim" _col(12) "attainable" _col(26) "visible" _col(40) "all"
            di as txt "{hline 44}"
            local nd = rowsof(sx_FACES)
            forvalues i = 1/`nd' {
                di as res _col(6) `=`i'-1' _col(14) %7.0f sx_FACES[`i', 1] ///
                   _col(27) %7.0f sx_FACES[`i', 2] _col(38) %7.0f sx_FACES[`i', 3]
            }
            di as txt "{hline 44}"
            di as txt "attainable = faces containing a projection attainable under"
            di as txt "  positive diagonal V; visible = faces of the facets visible"
            di as txt "  from the treated unit; all = faces of the donor hull."
        }
    }

    if "`allplacebos'" != "" {
        local pval = scalar(sx_pval)
        local rank = scalar(sx_rank)
        local n    = scalar(sx_n)
        di as txt _n "Permutation test (post/pre-MSPE ratio):"
        di as txt "  treated unit ranks " as res `rank' as txt " of " ///
           as res `n' as txt ",  permutation p-value = " as res %5.3f `pval'
    }

    if "`figure'" != "" {
        // variable labels for the default axis titles (the data are cleared below)
        local dlabel : variable label `depvar'
        if `"`dlabel'"' == "" local dlabel "`depvar'"
        local tlabel : variable label `time'
        if `"`tlabel'"' == "" local tlabel "`time'"
        preserve
        clear
        if "`allplacebos'" == "" {
            svmat sx_PATH, names(col)
            local gname = cond(strpos(`"`graphopts'"', "name(") == 0, "name(esc_path, replace)", "")
            // defaults, each replaced by the option of the same name when given
            if `"`title'"' == ""  local title `""Treated unit and synthetic control""'
            if `"`xtitle'"' == "" local xtitle `""`tlabel'""'
            if `"`ytitle'"' == "" local ytitle `""`dlabel'""'
            if `"`legend'"' == "" local legend `"order(1 "treated unit" 2 "synthetic control") ring(0) pos(2) col(1) region(lstyle(none))"'
            twoway (line treated year, lcolor(black) lwidth(medthick)) ///
                   (line synth year, lcolor(black) lpattern(dash)), ///
                   xline(`treatline', lpattern(dash) lwidth(medium) lcolor(gs6)) ///
                   legend(`legend') xtitle(`xtitle') ytitle(`ytitle') ///
                   title(`title') `gname' `graphopts'
        }
        else {
            svmat sx_GAP
            local nu = colsof(sx_GAP) - 1
            local tc = `=scalar(sx_trcol) + 1'
            // gap column k of sx_GAP is unit row k-1 of sx_INF (unit, pre-RMSE, post-RMSE, ratio)
            local trpre = sx_INF[`=`tc' - 1', 2]
            local cmd ""
            local nomit 0
            forvalues k = 2/`=`nu' + 1' {
                if `k' == `tc' continue
                if "`maxpre'" != "" {
                    if sx_INF[`=`k' - 1', 2] > `maxpre' * `trpre' {
                        local ++nomit
                        continue
                    }
                }
                local cmd `cmd' (line sx_GAP`k' sx_GAP1, lcolor(gs13) lwidth(vthin))
            }
            local cmd `cmd' (line sx_GAP`tc' sx_GAP1, lcolor(black) lwidth(medthick))
            local gname = cond(strpos(`"`graphopts'"', "name(") == 0, "name(esc_placebos, replace)", "")
            if `"`title'"' == ""  local title `""Placebo gaps (treated in black)""'
            if `"`xtitle'"' == "" local xtitle `""`tlabel'""'
            if `"`ytitle'"' == "" local ytitle `""gap: treated minus synthetic (`dlabel')""'
            if `"`legend'"' == "" local legend "off"
            twoway `cmd', yline(0, lcolor(gs10)) ///
                xline(`treatline', lpattern(dash) lwidth(medium) lcolor(gs6)) legend(`legend') ///
                xtitle(`xtitle') ytitle(`ytitle') ///
                title(`title') `gname' `graphopts'
            if "`maxpre'" != "" {
                di as txt "(placebo figure: `nomit' unit(s) with pre-treatment RMSE above `maxpre' times the treated unit's omitted; the test uses all units)"
            }
        }
        restore
    }

    ereturn post `b'
    ereturn scalar loss   = `loss'
    if "`slack'" != "" ereturn scalar slack = `slack'
    ereturn local  regime = "`regime'"
    if `exact' != -1 ereturn scalar exact = `exact'
    ereturn local  cmd    = "synthx"
    capture ereturn matrix balance = `bal'
    ereturn local  plabels = "`plabels'"
    ereturn local  vclass  = "`vclass'"
    if "`allplacebos'" != "" {
        ereturn scalar pvalue = `pval'
        ereturn scalar rank   = `rank'
        ereturn scalar n      = `n'
        ereturn matrix placebos = sx_INF
        ereturn matrix gaps     = sx_GAP
    }
    if "`vweights'" != "" {
        if scalar(sx_vdiag_ok) == 1 {
            tempname vd
            matrix `vd' = sx_VDIAG
            ereturn matrix Vdiag = `vd'
        }
        if scalar(sx_vpd_ok) == 1 {
            tempname vp
            matrix `vp' = sx_VPD
            ereturn matrix Vpd = `vp'
            ereturn scalar pd_mineig = scalar(sx_pd_mineig)
        }
    }
    if "`faces'" != "" {
        if scalar(sx_faces_ok) == 1 {
            tempname fc
            matrix `fc' = sx_FACES
            ereturn matrix faces = `fc'
        }
    }
    capture scalar drop sx_loss sx_slack sx_pval sx_rank sx_n sx_pre sx_post sx_err sx_one_err sx_exact
    capture scalar drop sx_vdiag_ok sx_vpd_ok sx_faces_ok sx_pd_mineig
    capture matrix drop sx_W sx_BAL sx_VDIAG sx_VPD sx_FACES sx_PATH sx_INF sx_GAP
end
