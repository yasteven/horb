import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from make_agn_catalog import build_catalog, canonical, concentration, main, validate_weights


class ActivityCatalogTests(unittest.TestCase):
    def test_absolute_mass_concentration_preserves_sign_diagnostics(self):
        def measure(values):
            return concentration([{"signed_mass_msun": x} for x in values])
        balanced = measure([1, -1])
        self.assertEqual(balanced["effective_components"], 2)
        self.assertEqual(balanced["negative_abs_mass_fraction"], 0.5)
        self.assertEqual(balanced["mass_cancellation_fraction"], 1)
        self.assertEqual(measure([3, 0])["effective_components"], 1)
        self.assertEqual(measure([0])["effective_components"], "")
        for key, value in measure([1, -2, 3]).items():
            self.assertAlmostEqual(value, measure([3e10, 1e10, -2e10])[key])

    def test_name_aliases_do_not_use_fuzzy_matching(self):
        self.assertEqual(canonical("UGC 02487"), canonical("ugc2487"))
        self.assertNotEqual(canonical("NGC 4388"), canonical("NGC4389"))

    def test_mismatched_weight_normalization_is_rejected(self):
        summary = dict(galaxy="test", k=1, net_mass_msun=100,
                       positive_mass_msun=100, negative_mass_msun=0)
        row = dict(state="1s", a0_star_kpc=1, weight=1, ref_mass_msun=100, signed_mass_msun=100)
        validate_weights(summary, [row])
        with self.assertRaisesRegex(ValueError, "normalization"):
            validate_weights(summary, [{**row, "weight": 2}])
        with self.assertRaisesRegex(ValueError, "summary"):
            validate_weights({**summary, "net_mass_msun": 90}, [row])

    def test_catalog_coverage_limits_and_jet_evidence(self):
        rows, weights, points, info = build_catalog()
        by_name = {r["galaxy"]: r for r in rows}
        self.assertEqual(len(rows), 175)
        self.assertEqual(len(weights), 696)
        self.assertEqual(len(points), 1472)
        matched = [r for r in rows if r["catalog_coverage"] == "matched"]
        self.assertEqual(len(matched), 51)
        self.assertEqual(sum(r["match_method"] == "SIMBAD_alias" for r in rows), 13)
        self.assertEqual(sum(r["optical_class"] in ("Seyfert", "LINER") for r in rows), 21)
        self.assertEqual(sum(r["radio_core_limit"] == "detection" for r in rows), 18)
        self.assertEqual(sum(r["radio_core_limit"] == "upper" for r in rows), 33)
        self.assertEqual(sum(r["xray_limit"] == "detection" for r in rows), 31)
        self.assertEqual(by_name["UGC02487"]["activity_name"], "NGC 1167")
        self.assertAlmostEqual(float(by_name["UGC02487"]["log_core_Lnu_1p5_W_Hz_published"]), 23.473909)
        self.assertEqual(by_name["UGC02487"]["radio_core_rank_among_AGN_detections"], 1)
        self.assertEqual(by_name["NGC6503"]["radio_core_rank_among_detections"], "")
        self.assertEqual(by_name["NGC4217"]["jet_evidence"], "radio_origin_unresolved")
        self.assertEqual(by_name["DDO154"]["catalog_coverage"], "no_catalog_match")
        self.assertEqual(by_name["DDO154"]["optical_class"], "unclassified")

    def test_export_contains_real_rotation_points_and_full_precision_weights(self):
        with tempfile.TemporaryDirectory() as directory:
            main(["--out-dir", directory])
            with (Path(directory) / "rotation_points.csv").open() as stream:
                points = list(csv.DictReader(stream))
            self.assertEqual(len(points), 1472)
            self.assertEqual(len(set(r["galaxy"] for r in points)), 51)
            ugc = [r for r in points if r["galaxy"] == "UGC02487"]
            self.assertEqual(len(ugc), 17)
            self.assertEqual(float(ugc[-1]["R_kpc"]), 80.38)
            self.assertTrue(all(float(r["eV_kms"]) > 0 for r in points))


if __name__ == "__main__":
    unittest.main()
