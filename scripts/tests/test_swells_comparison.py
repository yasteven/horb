"""Check independent candidate selection and matched scientific inputs."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_swells import DATA, compare, load_family, select_family


class SwellsComparisonTests(unittest.TestCase):
    def test_rotation_selection_uses_chi2_and_lensing_selection_skips_invalid_branches(self):
        def candidate(chi2, rms, theta, branches=1):
            return {"rotation_chi2": chi2, "rotation_rms_kms": rms,
                    "effective_einstein_radius_arcsec": theta,
                    "maximum_outward_crossings_per_ray": branches}
        pool = [candidate(1., 4., None, 2), candidate(2., 2., .84),
                candidate(3., 3., .89, 2), candidate(4., 5., .90)]
        rotation, lensing = select_family(pool, .89, "signed")
        self.assertEqual(rotation["rotation_chi2"], 1.)
        self.assertEqual(lensing["effective_einstein_radius_arcsec"], .90)
        self.assertFalse(rotation["selection_uses_lensing"])
        self.assertTrue(lensing["selection_uses_lensing"])
        with self.assertRaisesRegex(ValueError, "no single-branch"):
            select_family(pool[:1], .89, "signed")

    def test_reference_selections_are_distinct_and_preserve_provenance(self):
        signed = load_family(DATA/"reference/j2141", False)
        positive = load_family(DATA/"reference/j2141_nonnegative", True)
        before = copy.deepcopy((signed, positive))
        result = compare(signed, positive)
        self.assertEqual((signed, positive), before)
        expected = [("positive_only", "best_rotation", .8411486462),
                    ("positive_only", "closest_lensing", .8979233216),
                    ("signed", "best_rotation", .4279124155),
                    ("signed", "closest_lensing", .8428780921)]
        for row, (family, selection, theta) in zip(result["rows"], expected):
            self.assertEqual((row["family"], row["selection"]), (family, selection))
            self.assertAlmostEqual(row["effective_einstein_radius_arcsec"], theta, places=7)
            self.assertEqual(row["selection_uses_lensing"], selection == "closest_lensing")
        self.assertEqual(result["pools"]["signed"]["saved_candidates"], 21)
        self.assertEqual(result["pools"]["positive_only"]["saved_candidates"], 21)

    def test_comparison_rejects_mismatched_baryons_and_wrong_family(self):
        signed = load_family(DATA/"reference/j2141", False)
        positive = load_family(DATA/"reference/j2141_nonnegative", True)
        positive[0]["baryons"][0]["mass"] *= 2
        with self.assertRaisesRegex(ValueError, "different baryons"):
            compare(signed, positive)
        with self.assertRaisesRegex(ValueError, "sign constraint"):
            load_family(DATA/"reference/j2141", True)


if __name__ == "__main__":
    unittest.main()
