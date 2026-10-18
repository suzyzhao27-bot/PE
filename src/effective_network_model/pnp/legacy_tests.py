import json
import math
import time
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve
from solver import Geometry,solve,linear_conductance,assemble,bernoulli,F,CSCALE

ROOT=Path(__file__).resolve().parent
METRICS={}


class Tests(unittest.TestCase):
    def test_sg_equilibrium_and_derivative(self):
        v=np.array([-20.,-1.,-.00001,0.,.00001,1.,20.])
        b,db=bernoulli(v);bm,_=bernoulli(-v)
        self.assertLess(float(np.max(np.abs(b-bm*np.exp(-v)))),1e-12)
        epsilon=1e-5
        fd=(bernoulli(v+epsilon)[0]-bernoulli(v-epsilon)[0])/(2*epsilon)
        self.assertLess(float(np.max(np.abs(fd-db))),1e-7)

    def test_analytic_jacobian(self):
        g=Geometry(cells_per_fiber=4)
        r,y,e,m,h=solve(30,1,-.003,geometry=g)
        rng=np.random.default_rng(34014)
        y=y+rng.normal(0,.01,len(y));d=rng.normal(0,1,len(y))
        args=(m,y,30,1,-.003,1.33e-9,2.03e-9,298.15,78.3,True,0.)
        res,jac,sc,_=assemble(*args)
        epsilon=1e-6
        plus=assemble(m,y+epsilon*d,*args[2:],jacobian=False)[0]
        minus=assemble(m,y-epsilon*d,*args[2:],jacobian=False)[0]
        analytical=jac@d;finite=(plus-minus)/(2*epsilon)
        error=float(np.max(np.abs(analytical-finite))/np.max(np.abs(analytical)))
        self.assertLess(error,1e-6)
        METRICS["analytical_Jacobian_directional_relative_error"]=error

    def test_chemical_Robin_jacobian_and_exchange(self):
        g=Geometry(cells_per_fiber=4)
        chem={"site_density_nm2":.05,"K_Na_per_M":0.,"K_Cl_per_M":10000.}
        r,y,e,m,h=solve(30,1,geometry=g,surface_model=chem)
        rng=np.random.default_rng(4391)
        y=y+rng.normal(0,.01,len(y));d=rng.normal(0,1,len(y));eps=1e-6
        args=(m,y,30,1,0.,1.33e-9,2.03e-9,298.15,78.3,True,0.)
        res,jac,sc,_=assemble(*args,surface_model=chem)
        rp=assemble(m,y+eps*d,*args[2:],jacobian=False,surface_model=chem)[0]
        rm=assemble(m,y-eps*d,*args[2:],jacobian=False,surface_model=chem)[0]
        error=float(np.max(abs((rp-rm)/(2*eps)-jac@d))/np.max(abs(jac@d)))
        self.assertLess(error,1e-6)
        swapped=dict(chem,K_Na_per_M=10000.,K_Cl_per_M=0.)
        a=solve(30,1,geometry=g,surface_model=chem)[0]
        b=solve(30,1,geometry=g,surface_model=swapped,dp=2.03e-9,dm=1.33e-9)[0]
        self.assertAlmostEqual(a["common_salt_flux_mol_m2_s"],b["common_salt_flux_mol_m2_s"],places=10)
        self.assertAlmostEqual(a["right_minus_left_theoretical_aqueous_voltage_mv"],-b["right_minus_left_theoretical_aqueous_voltage_mv"],places=8)
        self.assertGreater(a["wall_sigma_min_c_m2"],-.05*1e18*1.602176634e-19)
        self.assertLess(a["wall_sigma_max_c_m2"],0)
        METRICS["chemical_Robin_Jacobian_directional_relative_error"]=error

    def test_exact_neutral_equal_diffusivity_laplace(self):
        g=Geometry(cells_per_fiber=8)
        r,y,e,m,h=solve(30,1,0,geometry=g,dp=1.6e-9,dm=1.6e-9)
        rhs=np.zeros(m.n);np.add.at(rhs,m.left,2*30);np.add.at(rhs,m.right,2*1)
        c=spsolve(m.L,rhs)
        reference=1.6e-9*CSCALE*float(sum(2*(c[m.right]-1)))/m.width_m
        self.assertLess(abs(r["common_salt_flux_mol_m2_s"]-reference)/reference,1e-8)
        self.assertLess(abs(r["right_minus_left_theoretical_aqueous_voltage_mv"]),1e-8)
        self.assertLess(float(np.max(abs(e["cp"]-c))),1e-8)
        self.assertLess(float(np.max(abs(e["cm"]-c))),1e-8)
        METRICS["neutral_geometry_Laplace_reference_flux_mol_m2_s"]=reference

    def test_charge_and_ion_exchange(self):
        g=Geometry(cells_per_fiber=8)
        neg=solve(30,1,-.00674949,geometry=g)[0]
        swapped=solve(30,1,.00674949,geometry=g,dp=2.03e-9,dm=1.33e-9)[0]
        self.assertAlmostEqual(neg["common_salt_flux_mol_m2_s"],swapped["common_salt_flux_mol_m2_s"],places=10)
        self.assertAlmostEqual(neg["right_minus_left_theoretical_aqueous_voltage_mv"],-swapped["right_minus_left_theoretical_aqueous_voltage_mv"],places=8)

    def test_staggered_multiple_fibers_conservation(self):
        g=Geometry(cells_per_fiber=4,fiber_rows=5)
        r,y,e,m,h=solve(30,1,-.00674949,geometry=g)
        self.assertEqual(int(np.sum(~m.mask)),5*g.cells_per_fiber**2)
        self.assertAlmostEqual(r["representative_network_span_nm"],239.4)
        self.assertLess(r["relative_net_current_error"],1e-7)
        self.assertLess(r["plus_conservation_relative"],1e-7)
        self.assertLess(r["minus_conservation_relative"],1e-7)
        self.assertLess(r["poisson_gauss_balance_scaled"],1e-8)

    def test_open_outer_current_root(self):
        g=Geometry(cells_per_fiber=4)
        monolithic=solve(30,1,-.003,geometry=g)[0]
        def current(voltage):
            return solve(30,1,-.003,geometry=g,open_circuit=False,right_voltage_v=voltage)[0]["right_ionic_current_a_m2"]
        root=brentq(current,-.06,.06,xtol=2e-12)
        voltage=monolithic["right_minus_left_theoretical_aqueous_voltage_mv"]/1000
        self.assertLess(abs(root-voltage),1e-9)
        METRICS["independent_outer_zero_current_root_voltage_difference_v"]=abs(root-voltage)

    def test_static_boltzmann_equilibrium(self):
        r,y,e,m,h=solve(1,1,-.01,geometry=Geometry(cells_per_fiber=8),open_circuit=False)
        self.assertLess(float(np.max(abs(e["cp"]*np.exp(e["u"])-1))),1e-8)
        self.assertLess(float(np.max(abs(e["cm"]*np.exp(-e["u"])-1))),1e-8)
        self.assertLess(abs(r["plus_flux_right_mol_m2_s"]),1e-10)
        self.assertLess(abs(r["minus_flux_right_mol_m2_s"]),1e-10)

    def test_linear_bias_halving(self):
        a,_=linear_conductance(1,-.01,geometry=Geometry(cells_per_fiber=8),bias_v=1e-4)
        b,_=linear_conductance(1,-.01,geometry=Geometry(cells_per_fiber=8),bias_v=5e-5)
        error=abs(a["projected_area_ionic_conductance_s_m2"]-b["projected_area_ionic_conductance_s_m2"])/b["projected_area_ionic_conductance_s_m2"]
        self.assertLess(error,1e-5)
        self.assertLess(abs(a["linear_current_cation_fraction_t_plus"]-b["linear_current_cation_fraction_t_plus"]),1e-5)
        METRICS["G_central_bias_100_to_50uV_relative_difference"]=error

    def test_field_and_flux_conservation(self):
        worst=0.
        for cl,sigma in ((30,-.01),(30,.01),(270,-.00674949)):
            r,y,e,m,h=solve(cl,1,sigma,geometry=Geometry(cells_per_fiber=8))
            self.assertLess(r["max_scaled_equation_residual"],1e-8)
            self.assertLess(r["poisson_gauss_balance_scaled"],1e-8)
            self.assertLess(r["relative_net_current_error"],1e-7)
            self.assertLess(r["plus_conservation_relative"],1e-7)
            self.assertLess(r["minus_conservation_relative"],1e-7)
            worst=max(worst,r["max_scaled_equation_residual"],r["poisson_gauss_balance_scaled"])
            zz=np.nonzero(m.mask)[0]
            axial=(zz[m.ej]-zz[m.ei])==1
            for z,c,d,expected in ((1,e["cp"],1.33e-9,r["plus_flux_right_mol_m2_s"]),(-1,e["cm"],2.03e-9,r["minus_flux_right_mol_m2_s"])):
                v=z*(e["u"][m.ej]-e["u"][m.ei])
                q=bernoulli(v)[0]*c[m.ei]-bernoulli(-v)[0]*c[m.ej]
                cuts=[d*CSCALE*sum(q[axial&(zz[m.ei]==j)])/m.width_m for j in range(m.ny-1)]
                error=max(abs(vv-expected) for vv in cuts)/abs(expected)
                self.assertLess(error,1e-7)
        METRICS["worst_scaled_Poisson_or_equation_residual_test_cases"]=worst


if __name__=="__main__":
    started=time.perf_counter()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    report=dict(status="PASS" if result.wasSuccessful() else "FAIL",tests_run=result.testsRun,
                failures=len(result.failures),errors=len(result.errors),metrics=METRICS,runtime_seconds=time.perf_counter()-started,
                scope="finite aqueous 2D representative local fiber layer; local u=0; boundary sigma given; not real PE/DI/network validation")
    (ROOT/"test_report.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    raise SystemExit(0 if result.wasSuccessful() else 1)
