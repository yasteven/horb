"""Independent force checks and regression checks for the published figures."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad

SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
from horb_highlights import G, MODELS, LABELS, dz2_v2, model_v2, fit_baseline, metrics, signed_density_xz, load_density
from make_readme_highlights import validate_weights
from horb_qi import velocity as qi_velocity


class HighlightPhysics(unittest.TestCase):
    def test_dz2_matches_independent_rust_multipole_goldens_and_scale_law(self):
        r = np.array([.5, 1, 5, 10, 25])
        expected = [73.869328853264207, 275.8974661279248, 5600.430683727943, 20447.645268698514, 16259.64979783057]
        np.testing.assert_allclose(dz2_v2(r, 1, 1e11), expected, rtol=1e-11)
        np.testing.assert_allclose(dz2_v2(r*3, 3, 2e11), np.array(expected)*2/3, rtol=1e-11)
        self.assertAlmostEqual(float(dz2_v2(1e6, 1, 1e11))/(G*1e11/1e6), 1, places=8)

    def test_classical_forces_match_direct_density_integrals(self):
        scale, mass, reference = 15., 2.5e11, 80.
        for model, alpha in [("nfw", .18), ("burkert", .18), ("einasto", .05), ("einasto", .18), ("einasto", 1.7)]:
            def integrand(x):
                if model == "nfw": return x/(1+x)**2
                if model == "burkert": return x*x/((1+x)*(1+x*x))
                return x*x*np.exp(-2/alpha*x**alpha)
            def enclosed(r):
                return quad(integrand, 0, r/scale, epsabs=1e-100, epsrel=2e-11, limit=300)[0]
            radii = np.array([.001, .01, .1, 1, 5, 25, 80, 120])
            expected = np.array([G*mass/r*enclosed(r)/enclosed(reference) for r in radii])
            actual = model_v2(model, radii, np.log10([mass, scale, alpha]), reference)
            np.testing.assert_allclose(actual, expected, rtol=2e-8, err_msg=model)

    def test_refitting_recovers_synthetic_curves(self):
        radii = np.geomspace(2, 80, 35)
        baryons = G*2e10*radii/(radii+2)**2
        for model in ("single", "nfw", "burkert", "einasto"):
            params = np.log10([3e11, 1.2 if model == "single" else 15, .3])
            target = np.sqrt(baryons+model_v2(model, radii, params, 80))
            result = fit_baseline(model, radii, target, np.ones(len(radii))*5, baryons, 80)
            self.assertLess(result["rms_kms"], 1e-6, model)
            self.assertFalse(result["bound_reached"])

    def test_density_is_signed_and_rejects_wrong_fit_or_duplicate_grid(self):
        term = {"n": 1, "l": 0, "m": 0, "a0_star_kpc": 2, "weight": -3., "ref_mass_msun": 1e11}
        x, z = np.meshgrid([-1., 0, 1.], [-1., 0, 1.])
        rho = signed_density_xz(x, z, [term])
        self.assertAlmostEqual(rho[1, 1]/(-3e11/(np.pi*2**3)), 1.)
        self.assertTrue(np.all(rho < 0))
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"density.csv"
            data = np.column_stack([x.ravel(), z.ravel(), rho.ravel()])
            def write():
                np.savetxt(path, data, delimiter=",", header="x_kpc,z_kpc,rho_signed_msun_per_kpc3", comments="")
            write()
            _, _, loaded = load_density(path, [term])
            np.testing.assert_allclose(loaded, rho)
            with self.assertRaisesRegex(ValueError, "does not match"):
                load_density(path, [{**term, "weight": 3.}])
            data[1] = data[0]
            write()
            with self.assertRaises(ValueError):
                load_density(path, [term])

    def test_published_scores_recompute_from_six_plotted_curves(self):
        # These small, tracked records let CI audit the figures without requiring
        # a new exhaustive search or any ignored reports/ directory contents.
        for name in ("mw", "m31"):
            prefix = ROOT/"data/highlights"/name
            report = json.loads(Path(str(prefix)+"_fit.json").read_text())
            comparison = json.loads(Path(str(prefix)+"_comparison.json").read_text())
            curves = np.genfromtxt(str(prefix)+"_curves.csv", delimiter=",", names=True)
            obs = np.array(report["fitted_observations"])
            indices = np.searchsorted(curves["r_kpc"], obs[:, 0])
            np.testing.assert_array_equal(curves["r_kpc"][indices], obs[:, 0])
            for model in MODELS:
                score = metrics(curves[f"v_{model}_kms"][indices], obs[:, 1], obs[:, 2])
                for key in ("rms_kms", "chi2"):
                    self.assertAlmostEqual(score[key], comparison["scores"][model][key], places=8)
            np.testing.assert_allclose(curves["v_superposition_kms"][indices]**2, obs[:, 3]+obs[:, 4])
            validate_weights(str(prefix)+"_weights.csv", report["components"])
            with self.assertRaises(ValueError):
                validate_weights(str(prefix)+"_weights.csv", report["components"][:-1])
            self.assertEqual(report["config"]["search"], "exhaustive")
            self.assertEqual(report["config"]["criterion"], "chi2")
            self.assertEqual(report["config"]["ridge"], 0)
            self.assertEqual(report["evaluated_subsets"], 31931)
            self.assertEqual(set(comparison["scores"]), set(MODELS))
            for model in ("nfw", "burkert", "einasto"):
                result = comparison["scores"][model]
                expected = np.sqrt(obs[:, 3]+model_v2(model, obs[:, 0], result["log10_parameters"], result["reference_radius_kpc"]))
                np.testing.assert_allclose(expected, curves[f"v_{model}_kms"][indices], rtol=1e-12)
            np.testing.assert_allclose(qi_velocity(obs[:, 0], obs[:, 3]), curves["v_qi_kms"][indices], rtol=1e-12)
            self.assertNotIn("v_single_kms", curves.dtype.names)
            self.assertEqual(comparison["scores"]["qi"]["parameter_count"], 0)
            self.assertFalse(comparison["scores"]["qi"]["fitted"])


if __name__ == "__main__":
    unittest.main()
