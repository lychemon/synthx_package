"""
Run the paper's known-answer solver cases against the package modules.
"""

import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "synthx_package" / "code" / "python"))

from synthx import simplex_solver, visible_facet_solver


def load_tests(loader, tests, pattern):
    path = ROOT / "synth-AS-replication" / "tests" / "test_solvers.py"
    if not path.is_file():
        def paper_cases_unavailable():
            raise unittest.SkipTest("The shared cases need the sibling synth-AS-replication directory.")
        tests.addTest(unittest.FunctionTestCase(paper_cases_unavailable))
        return tests
    spec = importlib.util.spec_from_file_location("paper_solver_regressions", path)
    cases = importlib.util.module_from_spec(spec)
    original_path = sys.path.copy()
    try:
        spec.loader.exec_module(cases)
    finally:
        sys.path[:] = original_path
    # Reuse the examples and expected answers with the package's implementation.
    cases.simplex_solver = simplex_solver
    cases.facets = visible_facet_solver
    tests.addTests(loader.loadTestsFromModule(cases))
    return tests
