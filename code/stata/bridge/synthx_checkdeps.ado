*! Check Python, solver dependencies, and a small Clarabel problem.

capture program drop synthx_checkdeps
program synthx_checkdeps, rclass
    version 16.1

    di as txt _n "{hline 60}"
    di as txt "synthx dependency check"
    di as txt "{hline 60}"

    di as txt "Stata version     : " as res c(stata_version)

    // Check that Stata can start its configured Python.
    capture python: pass
    if _rc {
        di as error _n "Python is NOT available to Stata."
        di as txt  "  Install Python 3.12 or later, then point Stata at it:"
        di as txt  `"      . python search"'
        di as txt  `"      . python set exec "C:/path/to/python.exe", permanently"'
        di as txt  "  Verify with:  . python query"
        return scalar ok = 0
        exit 198
    }
    di as txt "Python            : configured  (path via 'python query')"

    // Locate the bridge using the global path or adopath.
    if "$synthx_bridge" != "" {
        local pp "$synthx_bridge"
    }
    else {
        capture qui findfile synthx_preflight.py
        if _rc {
            di as error _n "Cannot find the synthx bridge folder (synthx_preflight.py)."
            di as txt  "  Make it visible to Stata in one of two ways:"
            di as txt  `"      . global synthx_bridge "<...>/code/stata/bridge""'
            di as txt  `"      . adopath ++ "<...>/code/stata/bridge""'
            return scalar ok = 0
            exit 601
        }
        local pp = subinstr("`r(fn)'", "synthx_preflight.py", "", .)
    }
    di as txt "bridge folder     : " as res "`pp'"
    di as txt "{hline 60}"

    // Remove a trailing backslash before passing the path to Python.
    while inlist(substr("`pp'", -1, 1), char(92), "/") {
        local pp = substr("`pp'", 1, length("`pp'") - 1)
    }
    // Load this installation's preflight even after another package was used.
    capture python: import os, importlib.util; _spec = importlib.util.spec_from_file_location("_synthx_preflight", os.path.join(r"`pp'", "synthx_preflight.py")); _preflight = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_preflight)
    if _rc {
        di as error _n "Could not load synthx_preflight.py from the bridge folder."
        di as txt  "  Confirm the bridge folder above is correct and readable."
        return scalar ok = 0
        exit 198
    }
    python: _ = _preflight.report()

    di as txt "{hline 60}"
    if scalar(sx_deps_ok) == 1 {
        di as result "All synthx dependency checks passed."
        return scalar ok = 1
    }
    else {
        di as error "Missing or incompatible dependencies (see report above)."
        di as txt  "  Use the Python executable shown above. From a terminal, run:"
        di as txt  `"    <that-python> -m pip install -r "<package>/code/python/requirements.txt""'
        di as txt  "  If the solver file is missing, check the synthx_bridge path."
        return scalar ok = 0
        exit 198
    }
    capture scalar drop sx_deps_ok
end
