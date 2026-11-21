import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from make_m31_data import DATA, G, observations, prepare, thin_disk_v2, baryonic_components
from horb_signed import baryons_function
from fit_signed_density import dataset_metadata


class AndromedaTests(unittest.TestCase):
    def test_published_rows_and_large_outer_errors_are_preserved(self):
        obs = observations(DATA/"corbelli2010_table1.csv")
        self.assertEqual(len(obs), 29)
        np.testing.assert_array_equal(obs[[0, 6, 14, 24, 28]], [
            [8.5, 232, 11], [14.5, 253.2, 5.6], [22.5, 223.2, 5.4],
            [32.5, 224.7, 32.1], [36.5, 240.1, 38.2]])
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            meta = prepare(tmp)
            normalized = np.genfromtxt(tmp/"rotation_curve.csv", delimiter=",", skip_header=1)
            np.testing.assert_array_equal(normalized, obs)
            self.assertIn("half the difference", meta["rotation_curve"]["errors"])
            self.assertIn("no separately resolved gas", " ".join(meta["baryons"]["notes"]))
            saved = (tmp/"dataset.json").read_bytes()
            prepare(tmp)
            self.assertEqual(saved, (tmp/"dataset.json").read_bytes())

    def test_thin_disk_matches_high_precision_reference_not_spherical_disk(self):
        # Independent 60-digit mpmath Bessel evaluation of Geehan Eq. 7.
        r = np.array([8.5, 14.5, 36.5])
        expected = [20900.9705907008252362, 21471.8095331505659071, 9661.5585563939026256]
        actual = thin_disk_v2(r, 7.2e10, 5.4)
        np.testing.assert_allclose(actual, expected, rtol=2e-13)
        spherical = G*7.2e10*(1-np.exp(-r/5.4)*(1+r/5.4))/r
        self.assertGreater(abs(actual[0]/spherical[0]-1), 0.1)
        model = json.loads((DATA/"geehan2006_model.json").read_text())
        parts = baryonic_components(r, model)
        np.testing.assert_allclose(parts[:, 1], actual)
        np.testing.assert_allclose(parts[:, 0], G*3.2e10*r/(r+0.61)**2)
        np.testing.assert_allclose(parts[:, 2], G*5.6e7/r)

    def test_provenance_rejects_a_different_target_or_baryon_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            prepare(tmp)
            args = [tmp/"dataset.json", tmp/"rotation_curve.csv", tmp/"baryons.csv"]
            meta = dataset_metadata(*args)
            self.assertEqual(meta["galaxy"], "M31 (Andromeda)")
            for path in args[1:]:
                before = path.read_bytes()
                path.write_bytes(before+b"\n")
                with self.assertRaisesRegex(ValueError, "does not match"):
                    dataset_metadata(*args)
                path.write_bytes(before)

    def test_force_table_preserves_signed_components_and_interpolates_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"baryons.csv"
            path.write_text("r_kpc,v2_disk_kms2,v2_gas_kms2,v2_baryons_kms2\n1,16,-12,4\n3,32,-16,16\n")
            evaluate = baryons_function(path)
            np.testing.assert_array_equal(evaluate([1, 2, 3]), [4, 10, 16])
            self.assertNotAlmostEqual(float(evaluate([2])[0]), ((2+4)/2)**2)
            with self.assertRaisesRegex(ValueError, "outside"):
                evaluate([0.9])
            with self.assertRaisesRegex(ValueError, "outside"):
                evaluate([3.1])
            path.write_text("r_kpc,v2_gas_kms2,v2_baryons_kms2\n1,-12,-12\n3,-16,-16\n")
            self.assertEqual(float(baryons_function(path)([2])[0]), -14)
            path.write_text("r_kpc,v2_disk_kms2,v2_baryons_kms2\n1,16,4\n3,32,16\n")
            with self.assertRaisesRegex(ValueError, "do not sum"):
                baryons_function(path)
            path.write_text("r_kpc,v2_baryons_kms2\n3,16\n1,4\n")
            with self.assertRaisesRegex(ValueError, "increasing"):
                baryons_function(path)

    def test_dense_baryon_table_covers_guard_domain_and_observations(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            prepare(tmp)
            evaluate = baryons_function(tmp/"baryons.csv")
            r = observations(DATA/"corbelli2010_table1.csv")[:, 0]
            model = json.loads((DATA/"geehan2006_model.json").read_text())
            np.testing.assert_allclose(evaluate(r), baryonic_components(r, model).sum(axis=1), rtol=1e-14)
            self.assertTrue(np.all(np.isfinite(evaluate([0.5, 60]))))
            # Between table nodes the force interpolation remains accurate.
            between = np.linspace(8.507, 36.493, 120)
            exact = baryonic_components(between, model).sum(axis=1)
            np.testing.assert_allclose(evaluate(between), exact, rtol=5e-6)


if __name__ == "__main__":
    unittest.main()
