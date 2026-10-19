* Helpers for the exhibit do-files. Load with do synthx_tools.do.

* synthx_series W newvar, unit(varname) time(varname) outcome(varname)
*   Build the synthetic outcome using W, whose columns are donor IDs.
*   The result is shared within each period and missing if a positive-weight
*   donor has a missing outcome.
capture program drop synthx_series
program synthx_series
    syntax anything(name=args id="W newvar"), Unit(varname) Time(varname) Outcome(varname)
    gettoken W newvar : args
    local newvar = trim("`newvar'")
    confirm matrix `W'
    confirm new variable `newvar'
    tempvar w wy miss
    quietly gen double `w' = 0
    local names : colnames `W'
    local J = colsof(`W')
    forvalues j = 1/`J' {
        local nm : word `j' of `names'
        quietly replace `w' = `W'[1, `j'] if `unit' == `nm'
    }
    quietly gen double `wy' = `w' * `outcome'
    quietly bysort `time': egen double `newvar' = total(`wy')
    quietly bysort `time': egen double `miss' = total((`w' > 1e-12) & missing(`outcome'))
    quietly replace `newvar' = . if `miss' > 0
end

* synthx_weights W, unit(varname) [names(varname) tol(#)]
*   Print weights above tol() (default 1e-6), largest first.
*   Use names() or the unit's value labels for donor names.
capture program drop synthx_weights
program synthx_weights
    syntax anything(name=W id="W"), Unit(varname) [NAMes(varname) TOL(real 1e-6)]
    confirm matrix `W'
    local cols : colnames `W'
    local J = colsof(`W')
    local ulab : value label `unit'
    local left ""
    forvalues j = 1/`J' {
        if abs(`W'[1, `j']) > `tol' local left `left' `j'
    }
    di as txt _n "Donor weights (nonzero), largest first:"
    di as txt "{hline 48}"
    local total 0
    while "`left'" != "" {
        local jmax ""
        local wmax -1
        foreach j of local left {
            if `W'[1, `j'] > `wmax' {
                local wmax = `W'[1, `j']
                local jmax `j'
            }
        }
        local left : list left - jmax
        local id : word `jmax' of `cols'
        local lab ""
        if "`names'" != "" {
            quietly levelsof `names' if `unit' == `id', local(lab) clean
        }
        if "`lab'" == "" & "`ulab'" != "" local lab : label `ulab' `id'
        if "`lab'" == "" | "`lab'" == "`id'" local lab "unit `id'"
        di as txt "  " as res %-30s "`lab'" _col(36) as res %8.4f `wmax'
        local total = `total' + `wmax'
    }
    di as txt "{hline 48}"
    di as txt "  " %-30s "sum" _col(36) as res %8.4f `total'
end

* synthx_unitid namevar "name" localname
*   Store the numeric panel ID for "name" in localname.
capture program drop synthx_unitid
program synthx_unitid
    args namevar name lname
    quietly levelsof `_dta[_TSpanel]' if `namevar' == "`name'", local(ids)
    if "`ids'" == "" {
        di as error "synthx_unitid: no unit named `name'"
        exit 498
    }
    c_local `lname' `ids'
end
