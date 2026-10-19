* Build the Stata panels from the CSV files. Run from code/.
* Input: ../data/smoking_data.csv (Abadie et al., 2010) and
*        ../data/german_reunification.csv (Abadie et al., 2015).
* Output: ../data/smoking.dta and ../data/germany.dta, with labeled numeric
* unit IDs and a declared panel structure.
set more off

* California / Proposition 99.
import delimited "../data/smoking_data.csv", clear varnames(1) case(preserve) ///
    stringcols(1) numericcols(2 3 4 5 6 7)
rename state statename
encode statename, gen(state)          // 1 = Alabama, ..., 3 = California, ...
replace year = round(year)
label variable cigsale   "Cigarette sales per capita (packs)"
label variable lnincome  "Log per-capita state income"
label variable beer      "Per-capita beer consumption"
label variable age15to24 "Share of population aged 15-24"
label variable retprice  "Retail cigarette price"
order state statename year
sort state year
xtset state year
compress
save "../data/smoking.dta", replace

* West Germany / reunification.
import delimited "../data/german_reunification.csv", clear varnames(1) ///
    case(preserve) stringcols(2) numericcols(1 3 4 5 6 7 8 9 10 11)
rename country countryname
encode countryname, gen(country)      // alphabetical; West Germany is last (17)
drop code
label variable gdp       "Per-capita GDP (current USD)"
label variable infrate   "Inflation rate"
label variable trade     "Trade openness"
label variable schooling "Schooling"
label variable invest60  "Investment rate (1960s average)"
label variable invest70  "Investment rate (1970s average)"
label variable invest80  "Investment rate (1980s average)"
label variable industry  "Industry share"
order country countryname year
sort country year
xtset country year
compress
save "../data/germany.dta", replace

di as result _n "00_prepare_data.do DONE: data/smoking.dta, data/germany.dta"
