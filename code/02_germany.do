* Reproduce West Germany / reunification (Section 7 and Appendix D).
* Run from code/ after 00_prepare_data.do creates ../data/germany.dta.
*
* Exhibits produced (output/):
*   table5_germany_weights.tex      Table 5: exact donor weights vs ADH (2015)
*   table6_predictor_weights_ger.tex Table 6: diagonal and PD predictor weights
*   tableD2_germany_pd.tex          Table D.2: diagonal vs PD donor weights
*   fig6_germany.{png,pdf}          Figure 6: West Germany and two synthetic controls
*   figD2_germany_pd.{png,pdf}      Figure D.2: diagonal vs PD classes
*   germany_results.txt            numerical summary
*
* Predictors use 1971-1980, except schooling (1970, 1975) and investment (1980).
* Validation uses GDP in 1981-1990, after the predictor windows; 16 donors.
set more off

adopath ++ "stata/bridge"
global synthx_bridge "`c(pwd)'/stata/bridge"
do synthx_tools.do
capture python clear

use "../data/germany.dta", clear
synthx_unitid countryname "West Germany" wg
di as txt "treated unit: West Germany = id `wg'"

* 1. Diagonal donor weights and predictor weights for both classes.
synthx gdp gdp trade infrate industry schooling(1970 1975) invest70(1980), ///
    trunit(`wg') trperiod(1990) prewindow(1981 1990) xperiod(1971 1980) vweights
matrix W_diag = e(b)
matrix VD     = e(Vdiag)
matrix VP     = e(Vpd)
local  loss_diag = e(loss)
assert e(regime) == "exterior"
assert e(exact) == 1

assert reldif(`loss_diag', 45803.57) < 1e-4      // paper: 45,804

* 2. Full positive definite class (Appendix D)
synthx gdp gdp trade infrate industry schooling(1970 1975) invest70(1980), ///
    trunit(`wg') trperiod(1990) prewindow(1981 1990) xperiod(1971 1980) vclass(pd)
matrix W_pd   = e(b)
local  loss_pd = e(loss)
assert reldif(`loss_pd', 20286.57) < 1e-4        // paper: 20,287

* 3. Published weights of Abadie, Diamond, and Hainmueller (2015), Table 1
matrix W_adh = W_diag * 0
foreach s in "Austria 0.42" "USA 0.22" "Japan 0.16" "Switzerland 0.11" "Netherlands 0.09" {
    gettoken nm val : s
    synthx_unitid countryname "`nm'" id
    local names : colnames W_adh
    local k : list posof "`id'" in names
    matrix W_adh[1, `k'] = `val'
}

* 4. Synthetic series, post-treatment gaps, figures
synthx_series W_diag synth_diag, unit(country) time(year) outcome(gdp)
synthx_series W_pd   synth_pd,   unit(country) time(year) outcome(gdp)
synthx_series W_adh  synth_adh,  unit(country) time(year) outcome(gdp)
preserve
    keep if country == `wg'
    keep year gdp synth_diag synth_pd synth_adh
    gen double gap_diag = gdp - synth_diag
    gen double gap_pd   = gdp - synth_pd
    gen double gap_adh  = gdp - synth_adh
    foreach v in diag pd adh {
        quietly summarize gap_`v' if inrange(year, 1991, 2003), meanonly
        local mgap_`v' = r(mean)
    }
    save "../output/germany_series.dta", replace

    * Match the paper's black lines, boxed panels, and inset legends.
    local gstyle graphregion(color(white) margin(small)) ///
        plotregion(lcolor(black) lwidth(thin) margin(medsmall)) ///
        ylabel(, angle(horizontal) nogrid) xlabel(1960(10)2000, nogrid) xtitle("year")
    local lsolid lcolor(black) lwidth(medthick) lpattern(solid)
    local ldash  lcolor(black) lwidth(medthick) lpattern(dash)
    local ldot   lcolor(black) lwidth(medthick) lpattern(dot)
    local vline  xline(1990, lpattern(dash) lwidth(medium) lcolor(gs6))
    local leg    ring(0) cols(1) region(lstyle(none)) size(medsmall)

    * Figure 6: West Germany, the exact synthetic control, and the ADH (2015) one
    twoway (line gdp year, `lsolid') (line synth_diag year, `ldash') ///
           (line synth_adh year, `ldot'), ///
        `vline' `gstyle' ytitle("per-capita GDP (current USD)") ///
        legend(order(1 "West Germany" 2 "synthetic (exact)" ///
                     3 "synthetic (ADH 2015 weights)") position(11) `leg') ///
        name(f6a, replace)
    twoway (line gap_diag year, `ldash') (line gap_adh year, `ldot'), ///
        yline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) `vline' `gstyle' ///
        ytitle("gap (treated - synthetic)") ///
        legend(order(1 "gap, exact" 2 "gap, ADH 2015 weights") position(7) `leg') ///
        name(f6b, replace)
    graph combine f6a f6b, rows(1) xsize(11) ysize(4.5) graphregion(color(white)) ///
        name(fig6, replace)
    graph export "../output/fig6_germany.png", replace width(2200)
    graph export "../output/fig6_germany.pdf", replace

    * Figure D.2: diagonal vs PD classes
    twoway (line gdp year, `lsolid') (line synth_diag year, `ldash') ///
           (line synth_pd year, `ldot'), ///
        `vline' `gstyle' ytitle("per-capita GDP (current USD)") ///
        legend(order(1 "West Germany" 2 "synthetic (diagonal class)" ///
                     3 "synthetic (full PD class)") position(11) `leg') ///
        name(fD2a, replace)
    twoway (line gap_diag year, `ldash') (line gap_pd year, `ldot'), ///
        yline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) `vline' `gstyle' ///
        ytitle("gap (treated - synthetic)") ///
        legend(order(1 "gap, diagonal class" 2 "gap, full PD class") position(7) `leg') ///
        name(fD2b, replace)
    graph combine fD2a fD2b, rows(1) xsize(11) ysize(4.5) graphregion(color(white)) ///
        name(figD2, replace)
    graph export "../output/figD2_germany_pd.png", replace width(2200)
    graph export "../output/figD2_germany_pd.pdf", replace
restore

* 5. Tables
local plab1 "Per-capita GDP (1971-80)"
local plab2 "Trade openness (1971-80)"
local plab3 "Inflation rate (1971-80)"
local plab4 "Industry share (1971-80)"
local plab5 "Schooling (1970, 1975)"
local plab6 "Investment rate (1970s)"

tempfile donors
preserve
    keep country countryname
    duplicates drop
    drop if country == `wg'
    gen double w_diag = 0
    gen double w_pd   = 0
    gen double w_adh  = 0
    local names    : colnames W_diag
    local names_pd : colnames W_pd
    local J = colsof(W_diag)
    forvalues j = 1/`J' {
        local nm : word `j' of `names'
        quietly replace w_diag = W_diag[1, `j'] if country == `nm'
        quietly replace w_adh  = W_adh[1, `j']  if country == `nm'
        local k : list posof "`nm'" in names_pd
        if `k' > 0 quietly replace w_pd = W_pd[1, `k'] if country == `nm'
    }
    save `donors'
restore

* Table 5: exact vs ADH (2015), sorted by the exact weight
preserve
    use `donors', clear
    keep if w_adh > 0 | w_diag > 1e-4
    gsort -w_diag -w_adh
    tempname fh
    file open `fh' using "../output/table5_germany_weights.tex", write replace
    file write `fh' "\begin{tabular}{lcc}" _n "\hline" _n
    file write `fh' "Donor & Exact (this paper) & Abadie et al. (2015) \\" _n "\hline" _n
    forvalues i = 1/`=_N' {
        file write `fh' (countryname[`i']) " & " %5.3f (w_diag[`i']) " & " %4.2f (w_adh[`i']) " \\" _n
    }
    file write `fh' "\hline" _n "\end{tabular}" _n
    file close `fh'
restore

* Table D.2: diagonal vs PD, with loss and mean gap
preserve
    use `donors', clear
    keep if w_diag > 1e-4 | w_pd > 1e-4
    gsort -w_diag -w_pd
    tempname fh
    file open `fh' using "../output/tableD2_germany_pd.tex", write replace
    file write `fh' "\begin{tabular}{lcc}" _n "\hline" _n
    file write `fh' "Donor & Diagonal class & Full positive definite class \\" _n "\hline" _n
    forvalues i = 1/`=_N' {
        file write `fh' (countryname[`i']) " & " %5.3f (w_diag[`i']) " & " %5.3f (w_pd[`i']) " \\" _n
    }
    file write `fh' "\hline" _n
    file write `fh' "Validation loss (1981-1990) & " %7.0fc (`loss_diag') " & " %7.0fc (`loss_pd') " \\" _n
    file write `fh' "Mean post-treatment gap, 1991-2003 (USD) & " %7.0fc (`mgap_diag') " & " %7.0fc (`mgap_pd') " \\" _n
    file write `fh' "\hline" _n "\end{tabular}" _n
    file close `fh'
restore

* Table 6: diagonal analytic center (column 2 of e(Vdiag)) and both PD centers.
tempname fh
file open `fh' using "../output/table6_predictor_weights_ger.tex", write replace
file write `fh' "\begin{tabular}{lccc}" _n "\hline" _n
file write `fh' " & & \multicolumn{2}{c}{Full positive definite} \\" _n "\cline{3-4}" _n
file write `fh' "Predictor & Diagonal & Log-det maximizer & Analytic center \\" _n "\hline" _n
forvalues i = 1/6 {
    file write `fh' "`plab`i'' & " %5.3f (VD[`i', 2]) " & " %5.3f (VP[`i', 1]) " & " %5.3f (VP[`i', 2]) " \\" _n
}
file write `fh' "\hline" _n "\end{tabular}" _n
file close `fh'

* 6. Summary (screen + output/germany_results.txt)
capture log close germany
log using "../output/germany_results.txt", replace text name(germany)
di as txt _n "West Germany: validation loss, diagonal class = " as res %9.1f `loss_diag' as txt "   (paper 45,804)"
di as txt    "West Germany: validation loss, full PD class  = " as res %9.1f `loss_pd'   as txt "   (paper 20,287)"
di as txt    "Mean post-treatment gap 1991-2003, exact      = " as res %8.0f `mgap_diag' as txt "   (paper -1,391)"
di as txt    "Mean post-treatment gap 1991-2003, ADH (2015) = " as res %8.0f `mgap_adh'  as txt "   (paper -1,707)"
di as txt    "Mean post-treatment gap 1991-2003, full PD    = " as res %8.0f `mgap_pd'   as txt "   (paper -1,582)"
matrix list W_diag, format(%6.4f) title("Exact donor weights, diagonal class")
matrix list W_pd,   format(%6.4f) title("Exact donor weights, full PD class")
matrix list VD, format(%6.4f) title("Diagonal rationalizers: log-det / center / min / max (Table 6, col. 1 = center)")
assert reldif(VD[1, 2], 0.443) < 2e-3   // Table 6: diagonal analytic-center GDP share
matrix list VP, format(%6.4f) title("PD rationalizers: log-det / center / min / max (Table 6, cols. 2-3)")
log close germany
* Remove machine-specific log-path metadata, including wrapped lines.
python: import re, io; p = r"../output/germany_results.txt"; t = io.open(p, encoding="latin-1").read(); t = re.sub(r"(?m)^[ \t]*log:[^\r\n]*(?:\r?\n[ \t]*>[^\r\n]*)*(?:\r?\n|$)", "", t); io.open(p, "w", encoding="latin-1", newline="").write(t)
di as result _n "02_germany.do DONE"
