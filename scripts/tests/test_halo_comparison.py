"""Check comparison consumers, shared scoring and archived CSV compatibility."""
import contextlib
import csv
import importlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from horb_common import has_einasto, einasto_label, score


class ScoreTests(unittest.TestCase):
    def test_observations_outside_model_are_excluded(self):
        observed = np.array([(1, 999, 1), (5, 100, 2), (10, 104, 2), (90, 999, 1)],
                            dtype=[("R_kpc", float), ("Vobs_kms", float), ("eV_kms", float)])
        row = score("test", [5, 10, 80], [100, 100, 100], observed, 0, 100)
        self.assertEqual(row[1:4], (5, 80, 2))
        self.assertAlmostEqual(row[4], np.sqrt(8))
        self.assertEqual(row[5:], (4, 2))

    def test_invalid_data_and_empty_overlap_are_rejected(self):
        observed = np.array([(10, 100, 0)],
                            dtype=[("R_kpc", float), ("Vobs_kms", float), ("eV_kms", float)])
        with self.assertRaises(ValueError):
            score("test", [5, 20], [100, 100], observed, 5, 20)
        observed["eV_kms"] = 1
        for radii in ([20, 5], [5, 5], [5, np.nan], [5]):
            with self.assertRaises(ValueError):
                score("test", radii, np.ones(len(radii)), observed, 5, 20)
        with self.assertRaises(ValueError):
            score("test", [20, 30], [100, 100], observed, 20, 30)


class ComparisonPrograms(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.r = np.array([0.05, 5, 10, 25, 50, 80])
        self.header = ["r_kpc", "v_horb_kms", "v_piso_kms", "v_nfw_kms", "v_burkert_kms",
                       "m_horb_msun", "m_piso_msun", "m_nfw_msun", "m_burkert_msun",
                       "v_einasto_kms", "m_einasto_msun", "einasto_r_minus2_kpc",
                       "einasto_alpha", "einasto_rho_minus2_msun_per_kpc3"]
        rows = [[r, 100, 110, 130, 140, 1e11, 1e11, 1e11, 1e11,
                 120, 1e11, 15, 0.18, 1.4e6] for r in self.r]
        self.dm = self.write("dm.csv", self.header, rows)
        self.old = self.write("old.csv", self.header[:9], [r[:9] for r in rows])
        self.baryons = self.write("baryons.csv", ["kind", "mass_msun", "scale_kpc"],
                                 [["bulge", 1e10, 1]])
        radii = [0.01, 5, 10, 25, 50, 80, 95]
        self.rc = self.write("rc.csv", ["R_kpc", "Vobs_kms", "eV_kms"],
                            [[r, np.sqrt(120**2 + self.baryons_v2(r)), 2] for r in radii])
        # Deliberately narrower than the spherical curves to exercise the overlap.
        self.cuda = self.write("cuda.csv", ["r_kpc", "v_horb_cuda_kms"],
                              [[r, 115] for r in [5, 10, 25, 50]])
        self.total = self.write("total.csv", ["r_kpc", "v_disk_kms", "v_bulge_kms"],
                               [[r, 0, np.sqrt(self.baryons_v2(r))] for r in self.r])

    @staticmethod
    def baryons_v2(r):
        return 4.3009e-6 * 1e10 * r / (r + 1)**2

    def write(self, name, header, rows):
        path = self.root / name
        with path.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
        return path

    def command(self, name, args):
        return subprocess.run([sys.executable, str(SCRIPTS / (name + ".py")), *map(str, args)],
                              check=True, capture_output=True, text=True)

    def test_metadata_and_archive_warning(self):
        new = np.genfromtxt(self.dm, delimiter=",", names=True)
        self.assertTrue(has_einasto(new))
        self.assertIn("0.18", einasto_label(new))
        self.assertIn("15", einasto_label(new))
        old = np.genfromtxt(self.old, delimiter=",", names=True)
        with contextlib.redirect_stderr(io.StringIO()) as warning:
            self.assertFalse(has_einasto(old))
        self.assertIn("regenerate", warning.getvalue())

    def test_all_models_use_same_baryons_and_observation_range(self):
        for cuda, count, maximum, n in [(False, 6, 80, 5), (True, 7, 50, 4)]:
            args = ["--dm-compare-csv", self.dm, "--rc-csv", self.rc,
                    "--baryons-csv", self.baryons, "--r-max", "95"]
            if cuda:
                args += ["--cuda-csv", self.cuda]
            result = self.command("score_dm_comparison", args)
            rows = list(csv.DictReader(io.StringIO(result.stdout)))
            self.assertEqual(len(rows), count)
            for row in rows:
                self.assertEqual(int(row["n"]), n)
                self.assertEqual(float(row["r_max"]), maximum)
            einasto = next(row for row in rows if row["model"] == "Einasto+baryons")
            self.assertEqual(float(einasto["chi2"]), 0)

    def test_legacy_score_and_summary_include_einasto(self):
        result = self.command("score_total_model_comparison_vs_mw", [
            "--total-csv", self.total, "--dm-compare-csv", self.dm, "--mw-csv", self.rc])
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(int(row["n"]) == 5 for row in rows))
        result = self.command("summarize_dm_comparison", [self.dm])
        self.assertIn("Einasto,", result.stdout)

    def test_every_comparison_plot_handles_new_and_archived_csv(self):
        for dm in [self.dm, self.old]:
            cases = {
                "plot_dm_comparison": [dm],
                "plot_total_model_comparison_vs_standard_rc": [
                    "--dm-compare-csv", dm, "--rc-csv", self.rc, "--baryons-csv", self.baryons],
                "plot_cuda_total_vs_mw": ["--cuda-csv", self.cuda,
                    "--dm-compare-csv", dm, "--rc-csv", self.rc, "--baryons-csv", self.baryons],
                "plot_total_model_comparison_vs_mw": [
                    "--total-csv", self.total, "--dm-compare-csv", dm, "--mw-csv", self.rc],
            }
            for name, args in cases.items():
                with self.subTest(script=name, archive=dm == self.old):
                    module = importlib.import_module(name)
                    out = self.root / (name + ".png")
                    argv = [name, *map(str, args), "-o", str(out)]
                    with patch.object(sys, "argv", argv), patch.object(plt, "close"), \
                            contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        module.main()
                    fig = plt.gcf()
                    labels = fig.axes[0].get_legend_handles_labels()[1]
                    self.assertEqual(any("Einasto" in label for label in labels), dm == self.dm)
                    self.assertGreater(out.stat().st_size, 1000)
                    plt.close(fig)


if __name__ == "__main__":
    unittest.main()
