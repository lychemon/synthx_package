"""
Known-answer cases for the public estimator and predictor weights.
"""

from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import warnings

import cvxpy as cp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code" / "python"))

from synthx import estimate
from synthx import interior_pd, predictor_weights, visible_facet_solver
from synthx.gpcheck import general_position_check


class EstimatorTests(unittest.TestCase):
    def test_fixed_metric_uses_its_off_diagonal_entries(self):
        X0 = np.array([[0.0, 1.0], [0.0, 1.0]])
        V = np.array([[4.0, 1.0], [1.0, 2.0]])
        result = estimate(np.array([1.0, 0.0]), X0, np.array([0.0]),
                          np.array([[0.0, 1.0]]), V=V)
        # Along (t, t), the predictor objective is minimized at t = 5/8.
        np.testing.assert_allclose(result["weights"], [3 / 8, 5 / 8], atol=1e-8)
        self.assertAlmostEqual(result["loss"], (5 / 8) ** 2, places=8)
        self.assertEqual(result["V_class"], "fixed")

    def test_fixed_metric_preserves_full_precision_weights(self):
        share = 0.123456789
        result = estimate(np.array([share]), np.array([[0.0, 1.0]]),
                          np.array([share]), np.array([[0.0, 1.0]]), V=[1.0])
        np.testing.assert_allclose(result["weights"], [1 - share, share], rtol=0, atol=1e-9)
        self.assertLess(result["loss"], 1e-16)

    def test_fixed_metric_is_invariant_to_shared_offsets_and_scalar_metric_scale(self):
        for offset, metric in ((1e10, 1.0), (0.0, 1e-24)):
            with self.subTest(offset=offset, metric=metric):
                result = estimate(np.array([offset + 0.25]),
                                  np.array([[offset, offset + 1.0]]),
                                  np.array([0.25]), np.array([[0.0, 1.0]]), V=[metric])
                np.testing.assert_allclose(result["weights"], [0.75, 0.25], rtol=0, atol=1e-8)
                self.assertLess(result["loss"], 1e-12)

    def test_interior_weights_match_the_treated_predictor(self):
        for metric in ("diagonal", "pd"):
            with self.subTest(metric=metric):
                result = estimate(np.array([0.25]), np.array([[0.0, 1.0, 2.0]]),
                                  np.array([0.5]), np.array([[0.0, 2.0, 5.0]]),
                                  V_class=metric)
                self.assertEqual(result["regime"], "interior")
                np.testing.assert_allclose(result["weights"], [0.75, 0.25, 0.0], atol=1e-8)
                self.assertLess(result["loss"], 1e-12)
                self.assertIn("cell", result["info"])

    def test_pd_reduces_a_collinear_hull_when_treated_is_in_its_span(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            result = estimate(np.array([3.0, 6.0]),
                              np.array([[0.0, 1.0, 2.0], [0.0, 2.0, 4.0]]),
                              np.array([4.0]), np.array([[0.0, 1.0, 2.0]]),
                              V_class="pd")
        np.testing.assert_allclose(result["weights"], [0.0, 0.0, 1.0], atol=1e-8)
        self.assertAlmostEqual(result["loss"], 4.0, places=8)
        self.assertEqual(result["info"]["reduced_rank"], 1)

    def test_interior_pd_search_finds_a_known_zero_loss_cell(self):
        X0 = np.array([[0.0, 2.0, 2.5, 0.0], [0.0, 0.0, 1.0, 2.0]])
        result = estimate(np.array([0.6, 0.2]), X0, np.array([0.6]),
                          np.array([[0.0, 1.0, 2.0, 3.0]]), V_class="pd")
        # The Euclidean cell (0, 1, 3) gives these barycentric weights and zero loss.
        np.testing.assert_allclose(result["weights"], [0.6, 0.3, 0.0, 0.1], atol=1e-8)
        self.assertEqual(result["regime"], "interior")
        self.assertLess(result["loss"], 1e-12)

    def test_interior_pd_solver_failure_does_not_return_a_margin(self):
        X0 = np.array([[0.0, 2.0, 2.5, 0.0], [0.0, 0.0, 1.0, 2.0]])
        with patch.object(cp.Problem, "solve", side_effect=cp.SolverError("solver failed")):
            with self.assertRaises((RuntimeError, cp.SolverError)):
                interior_pd.theta_pd((0, 1, 3), X0)

    def test_off_span_pd_metric_can_select_any_donor(self):
        X0 = np.array([[0.0, 1.0, 2.0], [0.0, 0.0, 0.0]])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            result = estimate(np.array([3.0, 1.0]), X0,
                              np.array([0.0]), np.array([[0.0, 1.0, 2.0]]),
                              V_class="pd")
        # V = [[1, -3], [-3, 10]] makes donor 0 the unique predictor fit.
        self.assertLess(result["loss"], 1e-9)
        np.testing.assert_allclose(result["weights"], [1.0, 0.0, 0.0], atol=1e-5)

    def test_diagonal_collinearity_reports_an_unsupported_hull(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            with self.assertRaises(ValueError):
                estimate(np.array([3.0, 6.0]),
                         np.array([[0.0, 1.0, 2.0], [0.0, 2.0, 4.0]]),
                         np.array([4.0]), np.array([[0.0, 1.0, 2.0]]))

    def test_exterior_results_keep_the_geometry_for_diagnostics(self):
        X0 = np.array([[0.0, 1.0, 3.0], [0.0, 2.0, 1.0]])
        X1 = np.array([0.0, 5.0])
        Z0, Z1 = np.array([[0.0, 1.0, 2.0]]), np.array([1.0])
        for solve in (visible_facet_solver.solve_exterior_hullfree,
                      visible_facet_solver.solve_exterior_pd):
            with self.subTest(solver=solve.__name__):
                loss, weights, info = solve(X0, X1, Z1, Z0)
                self.assertLess(loss, 1e-9)
                self.assertTrue(set(np.flatnonzero(weights > 1e-7)) <= set(info["face"]))
                self.assertTrue(info["visible_facets"])
        result = estimate(X1, X0, Z1, Z0)
        self.assertIn("general_position", result["info"])
        self.assertIn("exact", result["info"]["general_position"])

    def test_missing_geometry_does_not_certify_general_position(self):
        X0 = np.array([[0.0, 1.0, 3.0], [0.0, 2.0, 1.0]])
        X1 = np.array([0.0, 5.0])
        visible = [(tuple(face), normal, offset) for face, (normal, offset)
                   in visible_facet_solver.visible_facets(X0, X1).items()]
        for face, facets in ((None, visible), ((1,), None)):
            with self.subTest(face=face, facets_present=facets is not None):
                result = general_position_check(X0, X1, face, facets)
                self.assertFalse(result["exact"])
                self.assertTrue(result["messages"])


class PredictorWeightTests(unittest.TestCase):
    def setUp(self):
        self.X0 = np.array([[0.0, 1.0, 2.0, 0.0], [0.0, 0.0, 0.0, -1.0]])
        self.X1 = np.array([0.5, 1.0])
        self.weights = np.array([0.5, 0.5, 0.0, 0.0])

    def test_pd_center_omits_a_slack_that_is_always_zero(self):
        result = predictor_weights.pd_rationalizers(
            self.X1, self.X0, self.weights, standardized=True)
        # V12 = 0; one inactive slack vanishes and the other equals V22.
        np.testing.assert_allclose(result["Vstar"], np.diag([0.5, 0.5]), atol=1e-5)
        np.testing.assert_allclose(result["center"], np.diag([1 / 3, 2 / 3]), atol=1e-4)
        np.testing.assert_allclose(result["lo"], [0.0, 0.0], atol=1e-7)
        np.testing.assert_allclose(result["hi"], [1.0, 1.0], atol=1e-7)
        self.assertEqual(result["implicit_eq"], 1)

    def test_diagonal_center_and_ranges_keep_full_precision(self):
        result = predictor_weights.diag_rationalizers(
            self.X1, self.X0, self.weights, standardized=True)
        np.testing.assert_allclose(result["logdet"], [0.5, 0.5], atol=1e-5)
        np.testing.assert_allclose(result["center"], [1 / 3, 2 / 3], atol=1e-4)
        np.testing.assert_allclose(result["lo"], [0.0, 0.0], atol=1e-7)
        np.testing.assert_allclose(result["hi"], [1.0, 1.0], atol=1e-7)
        self.assertEqual(result["implicit_eq"], 1)
        self.assertGreater(np.max(np.abs(result["center"] - np.round(result["center"], 4))), 1e-6)

    def test_face_counts_include_closed_regions_at_a_boundary_vertex(self):
        X0 = np.array([[0.0, 1.0, 3.0], [0.0, 2.0, 1.0]])
        counts = predictor_weights.count_faces(np.array([1.0, 3.0]), X0,
                                               standardized=True)
        # Both visible edges have a feasible closed region containing only donor 1.
        np.testing.assert_array_equal(counts, [[1, 3, 3], [2, 2, 3]])


if __name__ == "__main__":
    unittest.main()
