"""
Check Stata's Python environment before loading the solver.

This file uses only the standard library until check() tests each dependency.
It can also be run directly from Python.
"""

from __future__ import annotations

import importlib
import os
import platform
import sys

MIN_PYTHON = (3, 12)
MIN_VERSION = {
    "numpy": (2, 5, 1),
    "scipy": (1, 18, 0),
    "cvxpy": (1, 9, 3),
    "clarabel": (0, 11, 1),
}


def _solver_dir():
    """
    Return the code/python directory next to this bridge.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.realpath(os.path.join(here, "..", "..", "python"))


def _version_tuple(version):
    """
    Read the leading numbers in a package version.
    """
    result = []
    for part in str(version).split(".")[:3]:
        digits = ""
        for char in part:
            if not char.isdigit():
                break
            digits += char
        result.append(int(digits) if digits else 0)
    return tuple(result)


def check():
    """
    Return (ok, report lines) for the dependencies, solver, and a test QP.
    """
    ok = True
    lines = ["Python executable : %s" % sys.executable,
             "Python version    : %s" % platform.python_version()]
    if sys.version_info < MIN_PYTHON:
        ok = False
        lines.append("  -> FAIL: need Python >= 3.12")

    modules = {}
    for package, minimum in MIN_VERSION.items():
        try:
            module = importlib.import_module(package)
        except Exception as error:
            ok = False
            lines.append("%-18s: import FAILED (%s)" % (package, error))
            continue
        modules[package] = module
        version = getattr(module, "__version__", "?")
        if _version_tuple(version) < minimum:
            ok = False
            lines.append("%-18s: %s (need >= %s)" %
                         (package, version, ".".join(map(str, minimum))))
        else:
            lines.append("%-18s: %s (ok)" % (package, version))

    solver_dir = _solver_dir()
    estimator = os.path.join(solver_dir, "synthx", "sc_estimator.py")
    if not os.path.isfile(estimator):
        ok = False
        lines.append("solver            : NOT FOUND at %s" % estimator)
    elif ok:
        try:
            sys.path[:] = [p for p in sys.path if os.path.realpath(p) != solver_dir]
            sys.path.insert(0, solver_dir)
            for name in list(sys.modules):
                if name == "synthx" or name.startswith("synthx."):
                    del sys.modules[name]
            importlib.invalidate_caches()
            engine = importlib.import_module("synthx.sc_estimator")
            if os.path.normcase(os.path.realpath(engine.__file__)) != os.path.normcase(estimator):
                raise ImportError("solver loaded from %s" % engine.__file__)
            lines.append("solver            : %s (imports ok)" % estimator)

            cp = modules["cvxpy"]
            weights = cp.Variable(2)
            problem = cp.Problem(cp.Minimize(cp.sum_squares(weights - [0.25, 0.75])),
                                 [weights >= 0, cp.sum(weights) == 1])
            problem.solve(solver="CLARABEL")
            np = modules["numpy"]
            if problem.status != cp.OPTIMAL or not np.allclose(weights.value, [0.25, 0.75], atol=1e-5):
                raise RuntimeError("Clarabel test QP did not return the expected solution")
            lines.append("Clarabel test QP  : solved (ok)")
        except Exception as error:
            ok = False
            lines.append("solver check      : FAILED (%s)" % error)
    else:
        lines.append("solver            : present; fix dependencies before testing")
    return ok, lines


def report():
    """
    Print the report and set sx_deps_ok when running inside Stata.
    """
    ok, lines = check()
    for line in lines:
        print(line)
    print("synthx dependencies: " + ("PASS" if ok else "FAIL"))
    try:
        from sfi import Scalar
    except ImportError:
        pass
    else:
        Scalar.setValue("sx_deps_ok", int(ok))
    return ok


if __name__ == "__main__":
    sys.exit(0 if report() else 1)
