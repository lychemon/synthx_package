* Template for your own panel. Run from code/; replace lines marked <<<.
* The example uses California data prepared by 00_prepare_data.do.
*
* Use one row per unit and period, a numeric unit ID, integer time, and
* numeric outcomes and predictors. Declare the panel with xtset.
* Plain predictors use xperiod(); var(periods) uses the specified periods.
* prewindow() selects the outcomes used to choose the donor weights.
* For separate training and validation, keep predictor periods before prewindow().
set more off

* 1. Load the command.
adopath ++ "stata/bridge"
global synthx_bridge "`c(pwd)'/stata/bridge"
do synthx_tools.do
synthx_checkdeps                       // check Python and the required packages
capture python clear

* 2. Load your panel.
use "../data/smoking.dta", clear       // <<< your data
xtset state year                       // <<< your unit id and time variable
* Use encode unitname, gen(unit) to turn string names into labeled numeric IDs.
synthx_unitid statename "California" treated   // <<< or:  local treated = <id>

* 3. Fit the diagonal class.
synthx cigsale lnincome age15to24 retprice ///        <<< outcome, plain predictors
       beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///  <<< special predictors
    trunit(`treated') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) /// <<< design
    vweights figure treatline(1988)      // <<< vertical line; defaults to trperiod()
* e(b)      donor weights (columns named by unit id)   e(loss)    validation loss
* e(regime) exterior / interior                        e(balance) K x 3 predictor balance
* e(Vdiag)  K x 4 diagonal predictor weights: log-det maximizer, analytic center, min, max
* e(Vpd)    K x 4 positive definite predictor weights: log-det maximizer, center, min, max
matrix list e(b)
matrix list e(Vdiag)

* Build the synthetic outcome path with the saved donor weights.
matrix W = e(b)
synthx_series W synth, unit(state) time(year) outcome(cigsale)   // <<< unit/time/outcome
list year cigsale synth if state == `treated' & inrange(year, 1985, 1995), noobs

* 4. Fit the full positive definite class.
synthx cigsale lnincome age15to24 retprice ///
       beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///
    trunit(`treated') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) vclass(pd)

* 5. Refit each unit as treated for the placebo test.
* postwindow() sets the post-treatment period used in the MSPE ratio.
local placebos 1                       // 0 skips the placebo test
if `placebos' {
    synthx cigsale lnincome age15to24 retprice ///
           beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), ///
        trunit(`treated') trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) ///
        allplacebos postwindow(1989 2000) figure treatline(1988) maxpre(1.5)   // <<< maxpre: hide poorly fit placebos in the figure
    di as txt "permutation p-value = " as res %5.3f e(pvalue) ///
       as txt "  (rank " as res e(rank) as txt " of " as res e(n) as txt ")"
    * e(placebos): unit / pre-RMSE / post-RMSE / ratio;  e(gaps): year x unit gaps
}

di as result _n "03_your_own_data.do DONE"
