"""QI units/redshift, and observing/scoring consistency in the SWELLS plot."""
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from horb_qi import C_M_S, THETA_M, KPC_M, A_QI_M_S2, velocity, metadata
from horb_highlights import metrics
from fit_swells import read_inputs, force_grid, total_force, ROTATION_MODELS
from horb_swells.models import baryon_v2, observing_matrix


class QIHighlightTests(unittest.TestCase):
    def test_acceleration_units_and_redshift_are_the_published_equations(self):
        r = np.array([.03, 8., 100.])
        b = A_QI_M_S2*r*KPC_M/1e6
        np.testing.assert_allclose(velocity(r, b)/np.sqrt(b), 2**.25, rtol=2e-15)
        for z in (0., .138, 2.):
            r_m, v_bar = 2.5e20, 153000.
            g_bar = v_bar**2/r_m
            expected = np.sqrt(r_m*np.sqrt(g_bar*(g_bar+2*C_M_S**2*(1+z)/THETA_M)))/1000
            self.assertAlmostEqual(float(velocity(r_m/KPC_M, (v_bar/1000)**2, z)), expected, places=11)
            self.assertEqual(metadata(z)["parameter_count"], 0)
        np.testing.assert_array_equal(velocity([1., 10.], 0.), [0., 0.])
        for r, b, z in [(0, 1, 0), (-1, 1, 0), (1, -1, 0), (1, np.nan, 0),
                        (1, 1, -1), (1, 1, np.inf)]:
            with self.assertRaises(ValueError):
                velocity(r, b, z)

    def test_j2141_prediction_is_projected_once_and_scores_only_six_outer_points(self):
        meta, obs = read_inputs()
        report = json.loads((ROOT/"data/swells/reference/j2141/rotation_fit.json").read_text())
        grid = force_grid(1600)
        expected_v = np.zeros_like(grid)
        expected_v[1:] = velocity(grid[1:], baryon_v2(grid[1:], report["baryons"]), meta["lens_redshift"])
        np.testing.assert_allclose(total_force("qi", grid, report), expected_v**2, rtol=1e-13)
        P = observing_matrix(obs["r_arcsec"], grid, meta, step_arcsec=.02)
        prediction = P@expected_v
        table = np.genfromtxt(ROOT/"data/highlights/j2141_predictions.csv", delimiter=",", names=True)
        np.testing.assert_allclose(prediction, table["v_qi_los_kms"], rtol=1e-11, atol=1e-9)
        self.assertNotIn("v_single_los_kms", table.dtype.names)
        saved = json.loads((ROOT/"data/highlights/j2141_comparison.json").read_text())
        self.assertEqual(saved["scores"]["qi"]["redshift"], .138)
        mask = obs["use_for_fit"].astype(bool)
        self.assertEqual(int(mask.sum()), 6)
        for name in ROTATION_MODELS:
            score = metrics(table[f"v_{name}_los_kms"][mask], obs["v_los_kms"][mask], obs["error_kms"][mask])
            for key in ("rms_kms", "chi2"):
                self.assertAlmostEqual(score[key], saved["scores"][name][key], places=7)


if __name__ == "__main__":
    unittest.main()
