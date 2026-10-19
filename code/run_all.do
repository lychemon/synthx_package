* Run the Stata examples and paper exhibits from code/.
* Each do-file can also be run separately.
set more off

* Check the command, Python dependencies, and data files.
adopath ++ "stata/bridge"
global synthx_bridge "`c(pwd)'/stata/bridge"
confirm file "stata/bridge/synthx.ado"
confirm file "stata/bridge/stata_backend.py"
confirm file "python/synthx/sc_estimator.py"
confirm file "../data/smoking_data.csv"
confirm file "../data/german_reunification.csv"
synthx_checkdeps                 // checks Python, NumPy, SciPy, CVXPY, and Clarabel
capture mkdir "../output"

local run_template 1             // 0 skips the data template and its placebo test

do 00_prepare_data.do            // CSV -> data/smoking.dta, data/germany.dta
do example_california.do         // worked example; writes example_*.png
do 01_california.do              // Tables 1-4, D.1; Figures 5, D.1
do 02_germany.do                 // Tables 5-6, D.2; Figures 6, D.2
if `run_template' do 03_your_own_data.do

di as result _n "ALL DONE -- see ../output/"
