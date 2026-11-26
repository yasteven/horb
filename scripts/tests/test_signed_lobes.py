"""Physical invariants, independent projection and zero-orbital regressions."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.special import kv

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from horb_lobes import density, orbital_v2, normalization, net_mass_fraction, mass_fields, G
from horb_swells.lensing import orbital_surface_modes, FourierLens
from horb_signed import Search, load_basis
from make_agn_catalog import concentration


class SignedLobeTests(unittest.TestCase):
    def test_absolute_mass_net_mass_and_nodal_signs(self):
        for n,ell in ((1,0),(2,0),(3,2),(5,2),(6,2)):
            from horb_lobes import intervals
            radial_integral=sum(quad(lambda r:r*r*abs(float(density(r,1,n,ell,1,1))),a,b,
                                     epsabs=1e-10)[0] for a,b in zip(intervals(n,ell)[:-1],intervals(n,ell)[1:]))
            angular=quad(lambda mu:abs(float(density(1,mu,n,ell,1,1)/density(1,1,n,ell,1,1))),-1,1,
                         points=[-1/np.sqrt(3),1/np.sqrt(3)],epsabs=1e-10)[0]*2*np.pi
            self.assertAlmostEqual(radial_integral*angular,1.,places=8)
        self.assertGreater(density(1,1,3,2,1,1),0)
        self.assertLess(density(1,0,3,2,1,1),0)
        self.assertLess(density(3,1,2,0,1,1),0)
        self.assertEqual(net_mass_fraction(3,2),0.)
        self.assertLess(net_mass_fraction(2,0),0.)

    def test_exponential_monopole_force_and_zero_response(self):
        r=np.array([.05,1,5,30,300.])
        expected=G/r*(1-np.exp(-r)*(1+r+r*r/2))
        np.testing.assert_allclose(orbital_v2(r,1,0,1,1),expected,rtol=1e-9)
        np.testing.assert_array_equal(orbital_v2(r,2,1,1,1),0.)
        np.testing.assert_allclose(orbital_v2(2*r,3,2,2,2),orbital_v2(r,3,2,1,1),rtol=1e-12)

    def test_projected_psi_matches_independent_inclined_line_integral(self):
        r=np.array([.3,2.,8.,30.])
        ci=.35;si=np.sqrt(1-ci*ci)
        for n,ell in ((1,0),(2,0),(3,2)):
            t={"n":n,"l":ell,"m":0,"a0_star_kpc":1.,"ref_mass_msun":1.,"density_model":"psi"}
            modes=orbital_surface_modes(r,t,ci,order=256)
            for phi in (0,.37,np.pi/2):
                predicted=sum(m*np.cos(k*phi) for m,k in zip(modes,(0,2,4)))
                expected=[]
                for rr in r:
                    def integrand(z):
                        radius=np.hypot(rr,z)
                        mu=(rr*np.sin(phi)*si+z*ci)/radius
                        return float(density(radius,mu,n,ell,1,1))
                    expected.append(quad(integrand,-np.inf,np.inf,epsabs=1e-11,epsrel=1e-10)[0])
                np.testing.assert_allclose(predicted,expected,rtol=2e-6,atol=1e-10)
            if ell==0 and n==1:
                np.testing.assert_allclose(modes[0],r*kv(1,r)/(4*np.pi),rtol=1e-7)
        # Signed projected density must survive into the Hessian (trace=2κ).
        grid=np.geomspace(1e-5,300,3000)
        modes=orbital_surface_modes(grid,t,ci)
        lens=FourierLens(grid,modes)
        field=lens.field(r,0.)
        expected=orbital_surface_modes(r,t,ci).sum(axis=0)
        np.testing.assert_allclose((field[2]+field[4])/2,expected,rtol=2e-4,atol=1e-10)

    def test_null_is_baryons_and_net_zero_is_not_zero_concentration(self):
        b=np.array([4.,9.,16.])
        A=np.ones((3,2))
        s=Search(A,b,np.sqrt(b),np.ones(3),A,b,ridge=0,criterion="chi2")
        best=s.run(max_components=0,budget=1)[0]
        self.assertEqual(best.subset,())
        self.assertEqual(best.chi2,0.)
        t={"n":3,"l":2,"weight":2.,"ref_mass_msun":1e9,"density_model":"psi"}
        fields=mass_fields(t)
        self.assertEqual(fields["signed_mass_msun"],0.)
        self.assertEqual(fields["absolute_component_mass_msun"],2e9)
        self.assertEqual(concentration([fields])["largest_abs_mass_fraction"],1.)
        self.assertEqual(concentration([fields])["negative_abs_mass_fraction"],.5)
        self.assertEqual(concentration([fields])["mass_cancellation_fraction"],1.)
        self.assertEqual(concentration([])["largest_abs_mass_fraction"],"")

    def test_rust_force_snapshots_match_independent_scipy_evaluation(self):
        for name in ("mw","m31","agc114905"):
            report=json.loads((ROOT/f"data/highlights/{name}_psi_fit.json").read_text())
            snap=np.array(report["fitted_observations"])
            expected=sum((t["weight"]*orbital_v2(snap[:,0],t["n"],t["l"],t["a0_star_kpc"],t["ref_mass_msun"])
                          for t in report["components"]),np.zeros(len(snap)))
            np.testing.assert_allclose(snap[:,4],expected,rtol=2e-10,atol=1e-7)

    def test_committed_psi_reports_do_not_reinterpret_old_weights(self):
        for name in ("mw","m31"):
            old=json.loads((ROOT/f"data/highlights/{name}_fit.json").read_text())
            new=json.loads((ROOT/f"data/highlights/{name}_psi_fit.json").read_text())
            self.assertEqual(old["model"],"signed_linear_orbital_density")
            self.assertEqual(new["density_model"],"psi")
            self.assertEqual(new["n_observations"],old["n_observations"])
            self.assertTrue(all(t["m"]==0 and t["l"]%2==0 for t in new["components"]))
            self.assertTrue(all(t["density_model"]=="psi" for t in new["components"]))


if __name__=="__main__":unittest.main()
