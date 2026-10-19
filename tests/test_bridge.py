"""
Data preparation and repeated calls through the Stata Python bridge.
"""

import csv
import importlib.util
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "synthx_package"
sys.path.insert(0, str(PACKAGE / "code" / "python"))

from synthx.predictor_weights import standardize


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_bridge():
    return load_module("synthx_test_bridge", PACKAGE / "code" / "stata" / "bridge" / "stata_backend.py")


def small_panel():
    share = 0.123456789
    return {
        "id": np.repeat([30.0, 10.0, 40.0, 20.0], 3),
        "year": np.tile([2000.0, 2001.0, 2002.0], 4),
        "x": np.repeat([1.0, share, 2.0, 0.0], 3),
        "y": np.array([2.0, 3.0, 4.0, 2 * share, 1 + 2 * share, 2 + 2 * share,
                       4.0, np.nan, 6.0, 0.0, 1.0, 2.0]),
    }


def fit_small_panel(bridge):
    return bridge.fit(small_panel(), "y", ["x"], [], "id", "year", 10,
                      2000, 2001, 2000, 2001, vclass="equal")


class BridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bridge = load_bridge()

    def test_incomplete_donors_are_dropped_without_changing_donor_order(self):
        result = fit_small_panel(self.bridge)
        np.testing.assert_array_equal(result["donors"], [20.0, 30.0])
        self.assertEqual(result["n_dropped"], 1)
        np.testing.assert_allclose(result["weights"], [0.876543211, 0.123456789], rtol=0, atol=1e-9)
        self.assertLess(result["loss"], 1e-16)

    def test_placebo_uses_the_same_fixed_metric_and_outcome_windows(self):
        data = small_panel()
        years = np.array([2000, 2001, 2002])
        paths = self.bridge.unit_series(data, "id", "year", "y", years)
        result = self.bridge.placebo_row(
            data, "y", ["x"], [], "id", "year", 20, 2000, 2001,
            2000, 2001, years, paths, 2002, 2002, "equal")
        np.testing.assert_array_equal(result["donors"], [10.0, 30.0])
        np.testing.assert_allclose(result["weights"], [1.0, 0.0], atol=1e-8)
        np.testing.assert_allclose(result["gap"], -2 * 0.123456789, atol=1e-8)
        self.assertAlmostEqual(result["ratio"], 1.0, places=8)

    def test_reloading_the_bridge_preserves_other_projects_modules(self):
        names = ("simplex_solver", "visible_facet_solver", "predictor_weights")
        sentinels = {name: types.ModuleType(name) for name in names}
        with patch.dict(sys.modules, sentinels):
            first = fit_small_panel(load_bridge())
            second = fit_small_panel(load_bridge())
            for name in names:
                self.assertIs(sys.modules[name], sentinels[name])
            for name in names:
                path = Path(sys.modules["synthx." + name].__file__).resolve()
                self.assertEqual(path.parent, PACKAGE / "code" / "python" / "synthx")
        np.testing.assert_allclose(first["weights"], second["weights"], atol=1e-12)
        self.assertAlmostEqual(first["loss"], second["loss"], places=12)

    def test_application_matrices_match_the_paper(self):
        source = ROOT / "synth-AS-replication" / "code" / "support" / "data_preparation.py"
        if not source.is_file():
            self.skipTest("The data comparison needs the sibling synth-AS-replication directory.")
        preparation = load_module(
            "paper_test_data", source)
        cases = [
            ("smoking_data.csv", "state", "California", "cigsale",
             ["lnincome", "age15to24", "retprice"],
             [("beer", list(range(1984, 1989))), ("cigsale", [1988]),
              ("cigsale", [1980]), ("cigsale", [1975])],
             (1970, 1988, 1980, 1988), preparation.build_california),
            ("german_reunification.csv", "country", "West Germany", "gdp",
             ["gdp", "trade", "infrate", "industry"],
             [("schooling", [1970, 1975]), ("invest70", [1980])],
             (1981, 1990, 1971, 1980), preparation.build_germany),
        ]
        for filename, unit, treated, outcome, reg, special, windows, build in cases:
            with self.subTest(application=treated):
                with open(PACKAGE / "data" / filename) as handle:
                    rows = list(csv.DictReader(handle))
                names = sorted({row[unit] for row in rows})
                ids = {name: i + 1 for i, name in enumerate(names)}
                variables = set(reg + [outcome] + [name for name, _ in special])
                data = {name: np.array([float(row[name]) if row[name] else np.nan
                                       for row in rows]) for name in variables}
                data["id"] = np.array([ids[row[unit]] for row in rows], dtype=float)
                data["year"] = np.array([float(row["year"]) for row in rows])
                result = self.bridge.build(data, outcome, reg, special, "id", "year",
                                           ids[treated], *windows)
                paper = build(ROOT / "synth-AS-replication" / "data" / filename)
                self.assertEqual([names[int(i) - 1] for i in result["donors"]], paper["donors"])
                for key in ("Z1", "Z0"):
                    np.testing.assert_array_equal(result[key], paper[key])
                np.testing.assert_array_equal(result["X1"], paper["X1_raw"])
                np.testing.assert_array_equal(result["X0"], paper["X0_raw"])
                X1, X0, _ = standardize(result["X1"], result["X0"])
                np.testing.assert_allclose(X1, paper["X1"], rtol=1e-13, atol=1e-12)
                np.testing.assert_allclose(X0, paper["X0"], rtol=1e-13, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
