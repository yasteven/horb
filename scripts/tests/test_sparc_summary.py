"""Numerical and sample-definition checks for the README overview figures."""
import json
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from make_sparc_summary import (MODELS, OUTPUT, aggregate, checked_inputs,
    classical_fits, comparison_rows, activity_rows, activity_statistics, fit_halo, halo_v2)
from horb_highlights import G
from horb_vms_mond import fit_vms, vms_v2, mond_velocity, MOND_A0_M_S2
from horb_qi import KPC_M


class SparcSummaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog, cls.curves, cls.fits, cls.fit_run = checked_inputs()

    def test_equal_galaxy_means_are_not_pooled_pointwise_errors(self):
        rows = [{"model": model, "galaxy": galaxy, "rms_kms": rms,
                 "n_points": n, "status": "ok"}
                for model in MODELS for galaxy, rms, n in [("A", 1., 100), ("B", 9., 1)]]
        for row in aggregate(rows):
            self.assertEqual(row["mean_rms_kms"], 5.)
            self.assertEqual(row["common_mean_rms_kms"], 5.)
        for row in rows:
            if row["galaxy"] == "B" and row["model"] == "qi":
                row.update(status="undefined", rms_kms=None)
        summary = {r["model"]:r for r in aggregate(rows)}
        self.assertEqual(summary["qi"]["galaxies"], 1)
        self.assertEqual(summary["nfw"]["galaxies"], 2)
        self.assertTrue(all(r["common_mean_rms_kms"] == 1. for r in summary.values()))

    def test_complete_inputs_and_no_silent_force_clipping(self):
        self.assertEqual(len(self.fits), 350)
        self.assertEqual(sum(len(c) for c in self.curves.values()), 3391)
        halos = classical_fits(self.curves)  # Recomputes cached scores from parameters.
        rows = comparison_rows(self.curves, self.fits, halos)
        self.assertEqual(len(rows), 1750)
        invalid = {(r["galaxy"],r["model"]) for r in rows if r["status"] != "ok"}
        self.assertEqual(invalid, {("UGC01281",m) for m in ("qi","mond","baryons")})
        self.assertTrue(all(r["rms_kms"] is None for r in rows if r["status"] != "ok"))
        self.assertEqual({r["common_galaxies"] for r in aggregate(rows)}, {174})

    def test_pseudo_isothermal_force_matches_density_integral(self):
        scale, mass, reference = 2.3, 1.7e10, 27.
        radii = np.array([1e-6, .1, 2., 17., reference])
        parameters = np.log10([mass, scale])
        integral = lambda upper: quad(lambda r: r*r/(1+(r/scale)**2),0,upper,epsabs=1e-30)[0]
        expected = np.array([G*mass/r*integral(r)/integral(reference) for r in radii])
        np.testing.assert_allclose(halo_v2("piso",radii,parameters,reference), expected, rtol=1e-10)
        for model in ("nfw","burkert","piso","einasto"):
            p = [*parameters, np.log10(.18)] if model=="einasto" else parameters
            self.assertAlmostEqual(float(halo_v2(model,[reference],p,reference)[0]),G*mass/reference,places=10)

    def test_halo_fit_recovers_a_known_rotation_curve(self):
        radius = np.linspace(.1, 15., 35)
        speed = np.sqrt(halo_v2("burkert",radius,np.log10([1.6e10,2.]),radius[-1]))
        curve = np.column_stack([radius,speed,np.linspace(1,3,len(radius)),np.zeros((len(radius),5))])
        fitted = fit_halo(("synthetic","burkert",curve))
        self.assertLess(fitted["rms_kms"],1e-6)
        self.assertAlmostEqual(fitted["parameters"]["scale_kpc"],2.,places=6)

    def test_activity_unknowns_limits_and_density_families_stay_distinct(self):
        rows = activity_rows(self.fits)
        stats = activity_statistics(rows)
        self.assertEqual(len(rows),102)
        for model in ("squared","psi"):
            s = stats[model]
            self.assertEqual((s["optical_AGN"],s["detections"],s["upper_limits"]),(21,12,9))
            self.assertEqual([s["groups"][g]["n"] for g in ("resolved_jets","other_optical_AGN","HII_comparison")],[7,14,29])
        by_key = {(r["galaxy"],r["density_model"]):r for r in rows}
        self.assertGreater(by_key["NGC2841","squared"]["concentration"],.87)
        self.assertLess(by_key["NGC2841","psi"]["concentration"],.42)

    def test_vms_fit_recovers_amplitude_and_mond_has_correct_limits(self):
        r = np.linspace(.5,18,30)
        b = np.linspace(-100,8000,len(r))
        target = np.sqrt(b+vms_v2(r,90000.,6.))
        fit = fit_vms(r,target,np.ones(len(r)),b)
        self.assertLess(fit["rms_kms"],1e-9)
        self.assertAlmostEqual(fit["amplitude_kms2"],90000.,places=6)
        transition_b = MOND_A0_M_S2*r*KPC_M/1e6
        np.testing.assert_allclose(mond_velocity(r,transition_b)/np.sqrt(transition_b),np.sqrt((1+np.sqrt(5))/2))
        self.assertAlmostEqual(float(mond_velocity(1e-10,40000.)),200.,places=8)
        with self.assertRaises(ValueError):
            mond_velocity(1.,-1.)


if __name__ == "__main__":
    unittest.main()
