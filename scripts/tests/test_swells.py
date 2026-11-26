"""Independent normalization, observing-geometry and lensing physics checks."""
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np
from scipy.integrate import quad, simpson

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from fit_swells import read_inputs, fingerprints, force_grid, total_force, MODELS
from horb_highlights import G, dz2_v2
from horb_swells.models import (distances, baryon_components, chameleon_v2, baryon_v2,
                               orbital_v2, dictionary, basis_matrix, observing_matrix)
from horb_swells.fitting import ProjectedSearch
from horb_swells.lensing import (orbital_surface_modes, spherical_surface, FourierLens,
                                chameleon_field, CombinedLens)


class SwellsPhysics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.meta, cls.obs = read_inputs()
        cls.distance = distances(cls.meta)
        cls.baryons = baryon_components(cls.meta)

    def test_exact_published_data_mask_and_distance_convention(self):
        self.assertEqual(len(self.obs), 10)
        np.testing.assert_array_equal(self.obs["v_los_kms"], [3.5,114.1,153.8,212.7,243.8,259.8,256.8,254.9,263.4,265.9])
        self.assertEqual(self.obs["use_for_fit"].sum(), 6)
        self.assertAlmostEqual(self.distance["kpc_per_arcsec"], 2.43773612544, places=9)
        self.assertNotIn("j2141_lensing.json", fingerprints())

    def test_orbital_force_matches_exact_1s_and_independent_rust_dz2_goldens(self):
        r = np.array([.5,1.,5.,10.,25.])
        expected = np.array([73.869328853264207,275.8974661279248,5600.430683727943,20447.645268698514,16259.64979783057])
        np.testing.assert_allclose(orbital_v2(r,3,2,1), expected, rtol=1e-11)
        enclosed = 1-np.exp(-2*r)*(1+2*r+2*r*r)
        np.testing.assert_allclose(orbital_v2(r,1,0,1), G*1e11*enclosed/r, rtol=1e-12)
        for t in dictionary(scales=(1.,)):
            actual = orbital_v2(np.array([1e5]),t["n"],t["l"],1.)
            np.testing.assert_allclose(actual, G*1e11/1e5, rtol=2e-7)

    def test_chameleon_density_mass_and_force_normalizations(self):
        for c in self.baryons:
            # Ellipsoidal volume element is 4 pi q3 r² dr.
            a, b = c["a"], c["b"]
            radial_mass = quad(lambda r: r*r*(b*b-a*a)/((r*r+a*a)*(r*r+b*b)), 0, np.inf)[0]
            coefficient = c["mass"]/(2*np.pi**2*c["q3"]*(b-a))
            self.assertAlmostEqual(4*np.pi*c["q3"]*coefficient*radial_mass/c["mass"],1.,places=9)
            far = 1e6
            self.assertAlmostEqual(float(chameleon_v2(np.array([far]),c)[0])/(G*c["mass"]/far),1.,places=5)

    def test_observing_operator_integrates_velocity_after_signed_force_sum(self):
        grid = force_grid(500)
        P = observing_matrix(self.obs["r_arcsec"],grid,self.meta)
        self.assertLess(np.max(np.abs(P[0])),1e-13)
        flat = P@np.full(len(grid),250.)
        self.assertTrue(np.all(flat[1:] < 250*np.sqrt(.96)))
        terms = dictionary(scales=(1.,))[:2]
        A, b = basis_matrix(grid,terms), baryon_v2(grid,self.baryons)
        expected = P@np.sqrt(b+.7*A[:,0]+.3*A[:,1])
        incorrect = np.sqrt(P@b+.7*(P@A[:,0])+.3*(P@A[:,1]))
        self.assertGreater(np.max(np.abs(expected[1:]-incorrect[1:])),1.)
        fine = observing_matrix(self.obs["r_arcsec"],grid,self.meta,step_arcsec=.02)
        self.assertLess(np.max(np.abs((P-fine)@np.sqrt(b))),.02)

    def test_projected_optimizer_recovers_synthetic_signed_weights(self):
        grid = force_grid(250)
        terms = dictionary(scales=(1.,))[:2]
        A, b = basis_matrix(grid,terms), baryon_v2(grid,self.baryons)
        mask = self.obs["use_for_fit"].astype(bool)
        P = observing_matrix(self.obs["r_arcsec"][mask],grid,self.meta)
        weights = np.array([1.2,-.2])
        observed = P@np.sqrt(b+A@weights)
        search = ProjectedSearch(A,b,P,observed,np.ones(len(observed)),max_weight=5)
        fit = search.fit((0,1))
        self.assertIsNotNone(fit)
        np.testing.assert_allclose(fit["weights"],weights,atol=3e-5)
        self.assertLess(fit["rms_kms"],1e-5)

    def test_projected_orbital_mass_and_signed_linearity(self):
        radius = np.geomspace(1e-5,10000,1800)
        for n,ell in ((1,0),(2,1),(3,2)):
            term = {"n":n,"l":ell,"a0_star_kpc":1.,"ref_mass_msun":1e11}
            modes = orbital_surface_modes(radius,term,.2)
            total = 2*np.pi*simpson(modes[0]*radius**2,x=np.log(radius))
            self.assertAlmostEqual(total/1e11,1.,places=7)
            lens = FourierLens(radius,modes/1e9)
            negative = FourierLens(radius,-2*modes/1e9)
            np.testing.assert_allclose(negative.field(2.,1.),-2*lens.field(2.,1.),rtol=1e-12)

    def test_deflection_derivatives_recover_surface_density(self):
        sigma = self.distance["sigma_critical_msun_kpc2"]
        for c in self.baryons:
            x,y,h = 1.7,.8,1e-4
            f = chameleon_field(x,y,c,sigma)
            dx = (chameleon_field(x+h,y,c,sigma)[:2]-chameleon_field(x-h,y,c,sigma)[:2])/(2*h)
            dy = (chameleon_field(x,y+h,c,sigma)[:2]-chameleon_field(x,y-h,c,sigma)[:2])/(2*h)
            np.testing.assert_allclose([dx[0],dx[1],dy[0],dy[1]],[f[2],f[3],f[3],f[4]],atol=2e-8)
            r = np.hypot(x,y/c["q"])
            density = c["mass"]/(2*np.pi*c["q"]*(c["b"]-c["a"]))*(1/np.hypot(r,c["a"])-1/np.hypot(r,c["b"]))
            self.assertAlmostEqual((f[2]+f[4])/2/(density/sigma),1.,places=11)

    def test_fourier_lens_matches_nfw_analytic_deflection(self):
        mass,scale,reference,sigma = 1e11,10.,13.,4e9
        radial = np.geomspace(1e-6,10000,4000)
        k = np.zeros((3,len(radial)))
        k[0] = spherical_surface(radial,"nfw",np.log10([mass,scale]),reference)/sigma
        lens = FourierLens(radial,k)
        x = np.array([.1,.5,1.,2.,10.])
        special = np.array([np.arccosh(1/t)/np.sqrt(1-t*t) if t<1 else
                            np.arccos(1/t)/np.sqrt(t*t-1) if t>1 else 1. for t in x])
        shape_ref = np.log1p(reference/scale)-(reference/scale)/(1+reference/scale)
        expected = mass/(np.pi*scale*shape_ref*sigma)*(np.log(x/2)+special)/x
        np.testing.assert_allclose(lens.field(x*scale,np.zeros_like(x))[0],expected,rtol=2e-6)

    def test_non_spherical_hessian_and_projection_convergence(self):
        radial = np.geomspace(1e-5,10000,3200)
        term = {"n":3,"l":2,"a0_star_kpc":1.,"ref_mass_msun":1e11}
        modes = orbital_surface_modes(radial,term,.2,192)
        finer = orbital_surface_modes(radial,term,.2,384)
        np.testing.assert_allclose(modes,finer,rtol=3e-4,atol=100.)
        lens = FourierLens(radial,modes/4e9)
        x,y,h = 1.7,.8,1e-4
        f = lens.field(x,y)
        dx = (lens.field(x+h,y)[:2]-lens.field(x-h,y)[:2])/(2*h)
        dy = (lens.field(x,y+h)[:2]-lens.field(x,y-h)[:2])/(2*h)
        np.testing.assert_allclose([dx[0],dx[1],dy[0],dy[1]],[f[2],f[3],f[3],f[4]],atol=2e-7)

    def test_spherical_gaussian_critical_curve_has_known_radius(self):
        radius = np.geomspace(1e-6,10000,4000)
        # kappa=2 exp(-r²/2); mean interior kappa =4(1-exp(-r²/2))/r².
        modes = np.zeros((3,len(radius)))
        modes[0] = 2*np.exp(-radius**2/2)
        lens = CombinedLens(FourierLens(radius,modes),[],1.)
        curve = lens.critical_curve(1.,angles=32)
        theta = curve["effective_einstein_radius_arcsec"]
        self.assertAlmostEqual(4*(1-np.exp(-theta*theta/2))/theta**2,1.,places=6)
        self.assertEqual(curve["maximum_outward_crossings_per_ray"],1)

    def test_multiple_critical_branches_are_not_reported_as_one_einstein_radius(self):
        radius = np.geomspace(1e-6,10000,100)
        lens = CombinedLens(FourierLens(radius,np.zeros((3,len(radius)))),[],1.)
        def eigenvalue(x,y):
            r = np.hypot(x,y)
            return (r-1)*(r-2)*(r-3)
        with patch.object(lens,"tangential_eigenvalue",side_effect=eigenvalue):
            curve = lens.critical_curve(1.,angles=16)
        self.assertEqual(curve["maximum_outward_crossings_per_ray"],2)
        self.assertIsNone(curve["effective_einstein_radius_arcsec"])
        self.assertFalse(curve["closed_star_shaped_outer_curve"])

    def test_recorded_fits_reconstruct_scores_and_keep_lensing_held_out(self):
        mask = self.obs["use_for_fit"].astype(bool)
        grid = force_grid()
        P = observing_matrix(self.obs["r_arcsec"],grid,self.meta)
        base = SCRIPTS.parent/"data/swells/reference"
        for folder in ("j2141","j2141_nonnegative","j2141_salpeter"):
            report = json.loads((base/folder/"rotation_fit.json").read_text())
            self.assertFalse(report["lensing_used_in_fit"])
            self.assertEqual(report["fit_inputs_sha256"],fingerprints())
            for name in (m for m in MODELS if m in report["models"]):
                prediction = P@np.sqrt(total_force(name,grid,report))
                residual = prediction[mask]-self.obs["v_los_kms"][mask]
                rms = np.sqrt(np.mean(residual**2))
                chi2 = np.sum((residual/self.obs["error_kms"][mask])**2)
                self.assertAlmostEqual(rms,report["models"][name]["rms_kms"],places=7)
                self.assertAlmostEqual(chi2,report["models"][name]["chi2"],places=6)


if __name__ == "__main__":
    unittest.main()
