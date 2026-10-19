* Reproduce California / Proposition 99 (Section 6 and Appendix D).
* Run from code/ after 00_prepare_data.do creates ../data/smoking.dta.
*
* Exhibits produced (output/):
*   table1_donor_weights.tex        Table 1: exact donor weights vs ADH (2010)
*   table2_predictor_weights.tex    Table 2: diagonal predictor weights
*   table3_faces.tex                Table 3: attainable / visible / all faces
*   table4_predictor_weights_pd.tex Table 4: PD log-det maximizer and center
*   tableD1_california_pd.tex       Table D.1: diagonal vs PD donor weights
*   fig5_california.{png,pdf}       Figure 5: California vs exact synthetic control
*   figD1_california_pd.{png,pdf}   Figure D.1: diagonal vs PD classes
*   california_results.txt         numerical summary
*
* Predictors use 1980-1988, except beer (1984-1988) and the listed outcome lags.
* Validation uses 1970-1988; treatment starts in 1989, with 38 donor states.
set more off
local faces 1          // 0 skips the exhaustive face count for Table 3

* Load the command and plotting helpers.
adopath ++ "stata/bridge"
global synthx_bridge "`c(pwd)'/stata/bridge"
do synthx_tools.do
capture python clear

use "../data/smoking.dta", clear
synthx_unitid statename "California" ca
di as txt "treated unit: California = id `ca'"

* 1. Diagonal donor weights, predictor weights, and face counts.
if `faces' local facesopt "faces"
synthx cigsale lnincome age15to24 retprice ///
       beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///
    trunit(`ca') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) ///
    vweights `facesopt'
matrix W_diag  = e(b)
matrix VD      = e(Vdiag)
matrix VP      = e(Vpd)
local  loss_diag = e(loss)
local  pd_mineig = e(pd_mineig)
if `faces' matrix FACES = e(faces)
assert e(regime) == "exterior"
assert e(exact) == 1

assert reldif(`loss_diag', 58.4566) < 1e-4       // paper: 58.46

* 2. Full positive definite class (Appendix D)
synthx cigsale lnincome age15to24 retprice ///
       beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///
    trunit(`ca') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) vclass(pd)
matrix W_pd   = e(b)
local  loss_pd = e(loss)
assert reldif(`loss_pd', 54.7035) < 1e-4         // paper: 54.70

* 3. Synthetic series, post-treatment gaps, figures
synthx_series W_diag synth_diag, unit(state) time(year) outcome(cigsale)
synthx_series W_pd   synth_pd,   unit(state) time(year) outcome(cigsale)
preserve
    keep if state == `ca'
    keep year cigsale synth_diag synth_pd
    gen double gap_diag = cigsale - synth_diag
    gen double gap_pd   = cigsale - synth_pd
    quietly summarize gap_diag if inrange(year, 1989, 2000), meanonly
    local mgap_diag = r(mean)
    quietly summarize gap_pd if inrange(year, 1989, 2000), meanonly
    local mgap_pd = r(mean)
    save "../output/california_series.dta", replace

    * Match the paper's black lines, boxed panels, and inset legends.
    local gstyle graphregion(color(white) margin(small)) ///
        plotregion(lcolor(black) lwidth(thin) margin(medsmall)) ///
        ylabel(, angle(horizontal) nogrid) xlabel(1970(5)2000, nogrid) xtitle("year")
    local lsolid lcolor(black) lwidth(medthick) lpattern(solid)
    local ldash  lcolor(black) lwidth(medthick) lpattern(dash)
    local ldot   lcolor(black) lwidth(medthick) lpattern(dot)
    local vline  xline(1988, lpattern(dash) lwidth(medium) lcolor(gs6))   // passage of Proposition 99
    local leg    ring(0) cols(1) region(lstyle(none)) size(medsmall)

    * Figure 5: California vs exact synthetic California, levels and gap
    twoway (line cigsale year, `lsolid') (line synth_diag year, `ldash'), ///
        `vline' `gstyle' ytitle("per-capita cigarette sales (packs)") ///
        legend(order(1 "California" 2 "synthetic California") position(1) `leg') ///
        name(f5a, replace)
    twoway (line gap_diag year, `lsolid'), ///
        yline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) `vline' `gstyle' ///
        ytitle("gap (treated - synthetic)") legend(off) name(f5b, replace)
    graph combine f5a f5b, rows(1) xsize(11) ysize(4.5) graphregion(color(white)) ///
        name(fig5, replace)
    graph export "../output/fig5_california.png", replace width(2200)
    graph export "../output/fig5_california.pdf", replace

    * Figure D.1: diagonal vs PD classes (levels and gaps)
    twoway (line cigsale year, `lsolid') (line synth_diag year, `ldash') ///
           (line synth_pd year, `ldot'), ///
        `vline' `gstyle' ytitle("per-capita cigarette sales (packs)") ///
        legend(order(1 "California" 2 "synthetic (diagonal class)" ///
                     3 "synthetic (full PD class)") position(7) `leg') ///
        name(fD1a, replace)
    twoway (line gap_diag year, `ldash') (line gap_pd year, `ldot'), ///
        yline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) `vline' `gstyle' ///
        ytitle("gap (treated - synthetic)") ///
        legend(order(1 "gap, diagonal class" 2 "gap, full PD class") position(1) `leg') ///
        name(fD1b, replace)
    graph combine fD1a fD1b, rows(1) xsize(11) ysize(4.5) graphregion(color(white)) ///
        name(figD1, replace)
    graph export "../output/figD1_california_pd.png", replace width(2200)
    graph export "../output/figD1_california_pd.pdf", replace
restore

* 4. Tables
* Predictor labels follow e(plabels).
local plab1 "Log per-capita state income"
local plab2 "Share of population aged 15-24"
local plab3 "Retail cigarette price"
local plab4 "Per-capita beer consumption"
local plab5 "Cigarette sales per capita (1988)"
local plab6 "Cigarette sales per capita (1980)"
local plab7 "Cigarette sales per capita (1975)"

* Collect donor weights for the tables.
tempfile donors
preserve
    keep state statename
    duplicates drop
    drop if state == `ca'
    gen double w_diag = 0
    gen double w_pd   = 0
    local names    : colnames W_diag
    local names_pd : colnames W_pd
    local J = colsof(W_diag)
    forvalues j = 1/`J' {
        local nm : word `j' of `names'
        quietly replace w_diag = W_diag[1, `j'] if state == `nm'
        local k : list posof "`nm'" in names_pd
        if `k' > 0 quietly replace w_pd = W_pd[1, `k'] if state == `nm'
    }
    * published weights of Abadie, Diamond, and Hainmueller (2010), Table 2
    gen double w_adh = 0
    foreach s in "Colorado 0.164" "Connecticut 0.069" "Montana 0.199" "Nevada 0.234" "Utah 0.334" {
        gettoken nm val : s
        quietly replace w_adh = `val' if statename == "`nm'"
    }
    save `donors'
restore

* Table 1: ADH (2010) vs exact, donors with positive weight in either column
preserve
    use `donors', clear
    keep if w_adh > 0 | w_diag > 1e-4
    sort statename
    tempname fh
    file open `fh' using "../output/table1_donor_weights.tex", write replace
    file write `fh' "\begin{tabular}{lcc}" _n "\hline" _n
    file write `fh' "State & Abadie et al. (2010) & Exact (this paper) \\" _n "\hline" _n
    forvalues i = 1/`=_N' {
        file write `fh' (statename[`i']) " & " %5.3f (w_adh[`i']) " & " %5.3f (w_diag[`i']) " \\" _n
    }
    file write `fh' "\hline" _n "\end{tabular}" _n
    file close `fh'
restore

* Table D.1: diagonal vs PD weights, losses, and mean gaps.
preserve
    use `donors', clear
    keep if w_diag > 1e-4 | w_pd > 1e-4
    gsort -w_diag
    tempname fh
    file open `fh' using "../output/tableD1_california_pd.tex", write replace
    file write `fh' "\begin{tabular}{lcc}" _n "\hline" _n
    file write `fh' "State & Diagonal class & Full positive definite class \\" _n "\hline" _n
    forvalues i = 1/`=_N' {
        file write `fh' (statename[`i']) " & " %5.3f (w_diag[`i']) " & " %5.3f (w_pd[`i']) " \\" _n
    }
    file write `fh' "\hline" _n
    file write `fh' "Validation loss (1970-1988) & " %5.2f (`loss_diag') " & " %5.2f (`loss_pd') " \\" _n
    file write `fh' "Mean post-treatment gap, 1989-2000 (packs) & " %5.1f (`mgap_diag') " & " %5.1f (`mgap_pd') " \\" _n
    file write `fh' "\hline" _n "\end{tabular}" _n
    file close `fh'
restore

* Table 2 uses the log-det weights in column 1 of e(Vdiag).
* Here they equal the analytic center because the retained predictors match exactly.
tempname fh
file open `fh' using "../output/table2_predictor_weights.tex", write replace
file write `fh' "\begin{tabular}{lc}" _n "\hline" _n "Predictor & Weight \\" _n "\hline" _n
forvalues i = 1/7 {
    file write `fh' "`plab`i'' & " %6.4f (VD[`i', 1]) " \\" _n
}
file write `fh' "\hline" _n "\end{tabular}" _n
file close `fh'

* Table 4: PD log-det maximizer and analytic center
file open `fh' using "../output/table4_predictor_weights_pd.tex", write replace
file write `fh' "\begin{tabular}{lcc}" _n "\hline" _n
file write `fh' "Predictor & Log-det maximizer & Analytic center \\" _n "\hline" _n
forvalues i = 1/7 {
    file write `fh' "`plab`i'' & " %5.3f (VP[`i', 1]) " & " %5.3f (VP[`i', 2]) " \\" _n
}
file write `fh' "\hline" _n "\end{tabular}" _n
file close `fh'

* Table 3: faces by dimension
if `faces' {
    file open `fh' using "../output/table3_faces.tex", write replace
    file write `fh' "\begin{tabular}{cccc}" _n "\hline" _n
    file write `fh' "Dimensions & \#Attainable Faces & \#Visible Faces & \#Faces \\" _n "\hline" _n
    forvalues i = 1/`=rowsof(FACES)' {
        file write `fh' "`=`i'-1'" " & " %5.0f (FACES[`i', 1]) " & " %5.0f (FACES[`i', 2]) " & " %5.0f (FACES[`i', 3]) " \\" _n
    }
    file write `fh' "\hline" _n "\end{tabular}" _n
    file close `fh'
}

* 5. Summary (screen + output/california_results.txt)
capture log close california
log using "../output/california_results.txt", replace text name(california)
di as txt _n "California: validation loss, diagonal class  = " as res %8.4f `loss_diag' as txt "   (paper 58.46)"
di as txt    "California: validation loss, full PD class   = " as res %8.4f `loss_pd'   as txt "   (paper 54.70)"
di as txt    "Mean post-treatment gap 1989-2000, diagonal  = " as res %6.1f `mgap_diag' as txt "   (paper -19.0)"
di as txt    "Mean post-treatment gap 1989-2000, full PD   = " as res %6.1f `mgap_pd'   as txt "   (paper -19.8)"
di as txt    "Smallest eigenvalue of the PD log-det maximizer = " as res %8.2e `pd_mineig' as txt "   (paper 3.7e-3)"
di as txt    "Upper bound on the standardized income share (PD class) = " as res %5.3f VP[1, 4] as txt "   (paper 0.038)"
assert reldif(VD[2, 1], 1/3) < 1e-3 & reldif(VD[3, 1], 1/3) < 1e-3 & reldif(VD[7, 1], 1/3) < 1e-3   // Table 2
matrix list W_diag, format(%6.4f) title("Exact donor weights, diagonal class")
matrix list W_pd,   format(%6.4f) title("Exact donor weights, full PD class")
matrix list VD, format(%6.4f) title("Diagonal rationalizers: log-det / center / min / max (Table 2 = log-det column)")
matrix list VP, format(%6.4f) title("PD rationalizers: log-det / center / min / max (Table 4)")
if `faces' matrix list FACES, title("Faces by dimension: attainable / visible / all (Table 3)")
log close california
* Remove machine-specific log-path metadata, including wrapped lines.
python: import re, io; p = r"../output/california_results.txt"; t = io.open(p, encoding="latin-1").read(); t = re.sub(r"(?m)^[ \t]*log:[^\r\n]*(?:\r?\n[ \t]*>[^\r\n]*)*(?:\r?\n|$)", "", t); io.open(p, "w", encoding="latin-1", newline="").write(t)
di as result _n "01_california.do DONE"
