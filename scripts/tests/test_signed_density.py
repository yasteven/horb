import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

SCRIPTS=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SCRIPTS))
from horb_signed import Search, load_basis


class SignedSearchTests(unittest.TestCase):
    def setUp(self):
        self.r=np.linspace(1,30,40)
        self.A=np.column_stack([1000*self.r/(self.r+1)**2,
                                1000*self.r/(self.r+10)**2,30+0*self.r])
        self.w=np.array([1.0,-10.0,0.0])
        self.b=np.full(len(self.r),5000.)
        self.v=np.sqrt(self.b+self.A@self.w)
        self.e=np.full(len(self.r),0.1)

    def search(self,nonnegative=False):
        return Search(self.A,self.b,self.v,self.e,self.A,self.b,ridge=0,nonnegative=nonnegative)

    def test_recovers_negative_density_weight_and_selects_component_count(self):
        self.assertTrue(np.any(self.v**2-self.b<0))
        best=self.search().run(3,"exhaustive",budget=10)[0]
        self.assertEqual(best.subset,(0,1))
        np.testing.assert_allclose(best.weights,self.w[:2],atol=1e-6)
        self.assertLess(best.chi2,1e-12)
        control=self.search(True).run(3,"exhaustive",budget=10)[0]
        self.assertGreater(control.chi2,100)

    def test_seeded_annealing_is_reproducible(self):
        a=self.search().run(2,"anneal",budget=7,seed=123)[0]
        b=self.search().run(2,"anneal",budget=7,seed=123)[0]
        self.assertEqual(a.subset,b.subset)
        np.testing.assert_array_equal(a.weights,b.weights)
        self.assertEqual(a.objective,b.objective)

    def test_total_force_constraints_include_radii_outside_fit(self):
        A=np.ones((4,1))*100
        s=Search(A,np.ones(4),np.ones(4)*10,np.ones(4),
                 np.array([[100.],[-100.]]),np.array([1.,1.]),ridge=0)
        c=s.fit((0,))
        self.assertIsNotNone(c)
        self.assertLessEqual(c.weights[0],0.010000001)
        self.assertGreaterEqual(c.min_total_v2,0)

    def test_null_model_and_search_budget(self):
        s=Search(self.A,self.b,np.sqrt(self.b),self.e,self.A,self.b)
        self.assertEqual(s.run(2,"beam")[0].subset,())
        with self.assertRaises(ValueError): self.search().run(3,"exhaustive",budget=4)
        with self.assertRaises(ValueError): self.search().run(2,"beam",budget=2)

    def test_outward_basis_force_is_not_squared_or_clipped(self):
        A=-np.arange(1,11,dtype=float)[:,None]*10
        b=np.full(10,1000.)
        v=np.sqrt(b+2*A[:,0])
        best=Search(A,b,v,np.full(10,0.01),A,b,ridge=0).run(1,"exhaustive",budget=2)[0]
        self.assertEqual(best.subset,(0,))
        self.assertAlmostEqual(best.weights[0],2.0,places=6)


class SignedCliTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.r=np.linspace(1,30,40)
        A=np.column_stack([1000*self.r/(self.r+1)**2,1000*self.r/(self.r+10)**2,30+0*self.r])
        headers=["schema_version","basis_id","state","n","l","m","a0_star_kpc","ref_mass_msun",
                 "gravity_backend","geometry","normalization","n_radial","n_mu","n_phi","softening_kpc","r_kpc","v2_signed_kms2"]
        rows=[]
        for i in range(3):
            for r,k in zip(self.r,A[:,i]):
                rows.append([1,i,f"{i+1}s",i+1,0,0,1,1e11,"cpu_multipole","azimuthal_mean","analytic_full_space",128,8,0,0,r,k])
        self.basis=self.write("basis.csv",headers,rows)
        self.baryons=self.write("baryons.csv",["kind","mass_msun","scale_kpc"],[["bulge",1e11,1]])
        b=4.3009e-6*1e11*self.r/(self.r+1)**2
        self.rc=self.write("rc.csv",["R_kpc","Vobs_kms","eV_kms"],
                           zip(self.r,np.sqrt(b+A@np.array([1.,-10.,0.])),np.full(len(self.r),0.02)))

    def write(self,name,header,rows):
        path=self.root/name
        with path.open("w",newline="") as f:
            w=csv.writer(f);w.writerow(header);w.writerows(rows)
        return path

    def command(self,script,args):
        return subprocess.run([sys.executable,str(SCRIPTS/script),*map(str,args)],capture_output=True,text=True)

    def test_fit_outputs_can_be_reconstructed_and_plotted(self):
        prefix=self.root/"fit"
        result=self.command("fit_signed_density.py",["--basis-csv",self.basis,"--rc-csv",self.rc,
            "--baryons-csv",self.baryons,"--out-prefix",prefix,"--r-min",1,"--r-max",30,
            "--search","exhaustive","--max-components",3,"--ridge",0,"--budget",10])
        self.assertEqual(result.returncode,0,result.stderr)
        report=json.loads(prefix.with_suffix(".json").read_text())
        self.assertEqual(report["k"],2)
        self.assertTrue(any(t["weight"]<0 for t in report["components"]))
        self.assertLess(report["chi2"],1e-10)
        basis=load_basis(self.basis)
        curve=np.genfromtxt(str(prefix)+"_curve.csv",delimiter=",",names=True)
        reconstructed=sum(t["weight"]*basis.at(curve["r_kpc"])[:,t["basis_id"]] for t in report["components"])
        np.testing.assert_allclose(reconstructed,curve["v2_dark_signed_kms2"],atol=1e-8)
        self.assertTrue(np.any(reconstructed<0))
        out=self.root/"rotation.png"
        plot=self.command("plot_signed_fit.py",["--fit-prefix",prefix,"-o",out])
        self.assertEqual(plot.returncode,0,plot.stderr)
        self.assertGreater(out.stat().st_size,1000)
        density=self.write("density.csv",["x_kpc","z_kpc","rho_signed_msun_per_kpc3"],
                           [(x,z,x*1e6) for z in [-1,0,1] for x in [-1,0,1]])
        out=self.root/"density.png"
        plot=self.command("plot_signed_density.py",[density,"-o",out])
        self.assertEqual(plot.returncode,0,plot.stderr)
        self.assertGreater(out.stat().st_size,1000)

    def test_rejects_speed_only_basis_and_extrapolation(self):
        bad=self.write("old.csv",["basis_id","r_kpc","v_cuda_kms"],[[0,1,2],[0,2,3]])
        with self.assertRaisesRegex(ValueError,"lost force signs"): load_basis(bad)
        basis=load_basis(self.basis)
        with self.assertRaises(ValueError): basis.at([0.5])
        with self.assertRaises(ValueError): basis.at([50])
        with self.assertRaises(ValueError): basis.at([np.nan])


if __name__=="__main__": unittest.main()
