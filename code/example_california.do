* Fit California's Proposition 99 example from Section 6. Run from code/.
* Input: ../data/smoking.dta; built by 00_prepare_data.do if missing.
* Output: example_*.png in ../output/.
* Predictors use 1980-1988, except beer (1984-1988) and the listed outcome lags.
* Validation uses 1970-1988; treatment starts in 1989, with 38 donor states.
set more off

* 1. Load the command and check its Python dependencies.
adopath ++ "stata/bridge"
global synthx_bridge "`c(pwd)'/stata/bridge"
do synthx_tools.do                       // helpers: synthx_series, synthx_weights, synthx_unitid
synthx_checkdeps

* 2. Load and declare the panel.
capture confirm file "../data/smoking.dta"
if _rc do 00_prepare_data.do
use "../data/smoking.dta", clear
xtset state year
synthx_unitid statename "California" ca   // numeric id of California (3)

* 3. Fit the diagonal class.
* Plain predictors use xperiod(); var(periods) uses the specified periods.
* prewindow() sets validation years. vweights adds predictor-weight summaries.
* figure plots outcomes; treatline(1988) marks the passage of Proposition 99.
synthx cigsale lnincome age15to24 retprice ///
       beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///
    trunit(`ca') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) ///
    vweights figure treatline(1988) ///
    title("California and its synthetic control") ytitle("per-capita cigarette sales (packs)") xtitle("year") ///
    legend(order(1 "California" 2 "synthetic California") ring(0) position(7) cols(1) region(lstyle(none))) ///
    graphopts(xlabel(1970(5)2000) ylabel(, angle(horizontal)))
graph export "../output/example_california_path.png", replace width(1800)

* Stored results
di as txt _n "regime: " as res e(regime) as txt "   validation loss: " as res %8.4f e(loss) ///
   as txt "   general position at the solution: " as res e(exact)
* Show positive donor weights by state name, largest first.
matrix W = e(b)
synthx_weights W, unit(state) names(statename)
* Show the log-det maximizer and analytic center. Columns 3-4 hold the ranges.
matrix VD = e(Vdiag)
matrix VP = e(Vpd)
matrix VW = VD[1..., 1..2], VP[1..., 1..2]
matrix colnames VW = diag:logdet diag:center pd:logdet pd:center
matrix list VW, format(%7.3f) title("predictor weights that rationalize the donor weights")

* 4. Build the synthetic outcome path and treatment gap.
synthx_series W synth, unit(state) time(year) outcome(cigsale)
gen double gap = cigsale - synth
quietly summarize gap if state == `ca' & inrange(year, 1989, 2000), meanonly
di as txt _n "mean gap 1989-2000 (packs per capita): " as res %6.1f r(mean)
list year cigsale synth gap if state == `ca' & inrange(year, 1986, 1994), noobs sep(0)

* 5. Fit the full positive definite class.
synthx cigsale lnincome age15to24 retprice ///
       beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///
    trunit(`ca') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) vclass(pd)
di as txt "validation loss, positive definite class: " as res %8.4f e(loss)

* 6. Refit each state as treated for the placebo test.
* e(pvalue) is California's rank by post/pre MSPE ratio divided by 39.
* postwindow() sets the post-treatment years. maxpre(1.5) hides states whose
* pre-treatment RMSE exceeds 1.5 times California's; they remain in the test.
synthx cigsale lnincome age15to24 retprice ///
       beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///
    trunit(`ca') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) ///
    allplacebos postwindow(1989 2000) figure treatline(1988) maxpre(1.5) ///
    title("Placebo gaps: California in black") ytitle("gap in per-capita cigarette sales (packs)") xtitle("year")
graph export "../output/example_california_placebos.png", replace width(1800)
di as txt _n "permutation p-value = " as res %5.3f e(pvalue) ///
   as txt "  (rank " as res e(rank) as txt " of " as res e(n) as txt ")"

di as result _n "example_california.do DONE"
