"""Targeted verification of prescribed common-ion energy in the v3 SG solver."""
import json
import math
import time
import unittest
from pathlib import Path
import numpy as np
import solver_v3_baseline as baseline
from solver import Geometry, CommonIonEnergy, solve, assemble, bernoulli

ROOT = Path(__file__).resolve().parent
METRICS = {}


def cut_flux_errors(result, extra, mesh):
    axial_row = np.nonzero(mesh.mask)[0]
    axial = axial_row[mesh.ej]-axial_row[mesh.ei] == 1
    errors = []
    for z,c,d,key in ((1,extra['cp'],1.33e-9,'plus_flux_right_mol_m2_s'),
                      (-1,extra['cm'],2.03e-9,'minus_flux_right_mol_m2_s')):
        v = z*(extra['u'][mesh.ej]-extra['u'][mesh.ei])+extra['common_energy_kbt'][mesh.ej]-extra['common_energy_kbt'][mesh.ei]
        flux = bernoulli(v)[0]*c[mesh.ei]-bernoulli(-v)[0]*c[mesh.ej]
        target = result[key]
        cuts = [d*sum(flux[axial & (axial_row[mesh.ei] == j)])/mesh.width_m for j in range(mesh.ny-1)]
        errors.append(float(max(abs(q-target) for q in cuts)/max(abs(target),1e-30)))
    return errors


class EnergyTests(unittest.TestCase):
    def test_field_definition_and_domains(self):
        mesh = Geometry(cells_per_fiber=8,fiber_rows=5).build()
        model = CommonIonEnergy('two_sided_exponential',3.,20.)
        w,metadata = model.build(mesh)
        x = (np.nonzero(mesh.mask)[0]+.5)*mesh.h
        xl,xr = metadata['resolved_x_left_nm'],metadata['resolved_x_right_nm']
        inside = (x>=xl)&(x<xr)
        self.assertAlmostEqual(xr-xl,239.4)
        self.assertTrue(np.all(w[~inside] == 0))
        expected = 3*(np.exp(-(x[inside]-xl)/20)+np.exp(-(xr-x[inside])/20))
        self.assertLess(float(max(abs(w[inside]-expected))),1e-14)
        self.assertEqual(len(w),mesh.n)
        with self.assertRaises(ValueError): CommonIonEnergy('two_sided_exponential',1.,0.).build(mesh)

    def test_zero_energy_regression_against_preserved_v3(self):
        maximum = 0.
        for rows,sigma,c in ((1,0.,30.),(1,-.006,270.),(5,-.006,30.)):
            kwargs = dict(c_left=c,c_right=1.,sigma=sigma)
            old = baseline.solve(**kwargs,geometry=baseline.Geometry(cells_per_fiber=4,fiber_rows=rows))
            for kind in ('rectangle','two_sided_exponential'):
                new = solve(**kwargs,geometry=Geometry(cells_per_fiber=4,fiber_rows=rows),ion_energy_model=CommonIonEnergy(kind,0.,20.))
                difference = float(max(abs(new[1]-old[1])))
                maximum = max(maximum,difference)
                self.assertLess(difference,1e-12)
                self.assertLess(abs(new[0]['common_salt_flux_mol_m2_s']-old[0]['common_salt_flux_mol_m2_s']),1e-12)
        METRICS['zero_energy_max_state_difference_from_v3'] = maximum

    def test_boltzmann_detailed_balance_including_jumps(self):
        records = []
        for kind in ('rectangle','two_sided_exponential'):
            for sigma in (0.,-.003):
                r,y,e,m,h = solve(1.,1.,sigma,geometry=Geometry(cells_per_fiber=8,fiber_rows=3),
                                  open_circuit=False,tol=1e-12,ion_energy_model=CommonIonEnergy(kind,4.,20.))
                w = e['common_energy_kbt']
                chemical_error = max(float(max(abs(e['cp']*np.exp(e['u']+w)-1))),
                                     float(max(abs(e['cm']*np.exp(-e['u']+w)-1))))
                self.assertLess(chemical_error,1e-8)
                local_flux = []
                for z,c in ((1,e['cp']),(-1,e['cm'])):
                    v=z*(e['u'][m.ej]-e['u'][m.ei])+w[m.ej]-w[m.ei]
                    local_flux.append(float(max(abs(bernoulli(v)[0]*c[m.ei]-bernoulli(-v)[0]*c[m.ej]))))
                self.assertLess(max(local_flux),1e-9)
                self.assertLess(abs(r['plus_flux_right_mol_m2_s']),1e-9)
                self.assertLess(abs(r['minus_flux_right_mol_m2_s']),1e-9)
                records.append(dict(kind=kind,sigma=sigma,boltzmann_relative_error=chemical_error,
                                    max_dimensionless_internal_face_flux=max(local_flux),residual=r['max_scaled_equation_residual']))
        METRICS['detailed_balance'] = records

    def test_analytic_jacobian_with_nonzero_energy(self):
        errors = []
        for kind in ('rectangle','two_sided_exponential'):
            r,y,e,m,h = solve(30.,1.,-.003,geometry=Geometry(cells_per_fiber=4,fiber_rows=3),
                              ion_energy_model=CommonIonEnergy(kind,2.,20.))
            rng=np.random.default_rng(1017);y=y+rng.normal(0,.01,len(y));d=rng.normal(0,1,len(y))
            args=(m,y,30.,1.,-.003,1.33e-9,2.03e-9,298.15,78.3,True,0.)
            kw=dict(common_energy=e['common_energy_kbt'])
            res,jac,sc,ex=assemble(*args,**kw);delta=1e-6
            rp=assemble(m,y+delta*d,*args[2:],jacobian=False,**kw)[0]
            rm=assemble(m,y-delta*d,*args[2:],jacobian=False,**kw)[0]
            error=float(max(abs((rp-rm)/(2*delta)-jac@d))/max(abs(jac@d)))
            self.assertLess(error,1e-6);errors.append(error)
        METRICS['energy_jacobian_directional_relative_errors'] = errors

    def test_positive_barrier_suppression_and_open_circuit_conservation(self):
        records=[];g=Geometry(cells_per_fiber=8,fiber_rows=5)
        for sigma in (0.,-.005):
            zero=solve(30.,1.,sigma,geometry=g)[0]
            for kind in ('rectangle','two_sided_exponential'):
                previous=zero['common_salt_flux_mol_m2_s']
                for amplitude in (1.,3.):
                    r,y,e,m,h=solve(30.,1.,sigma,geometry=g,ion_energy_model=CommonIonEnergy(kind,amplitude,26.6))
                    self.assertLess(r['common_salt_flux_mol_m2_s'],previous)
                    previous=r['common_salt_flux_mol_m2_s']
                    for key in ('relative_net_current_error','plus_conservation_relative','minus_conservation_relative',
                                'poisson_gauss_balance_scaled','max_scaled_equation_residual'):
                        self.assertLess(r[key],1e-7)
                    cut_error=max(cut_flux_errors(r,e,m));self.assertLess(cut_error,1e-7)
                    records.append(dict(kind=kind,amplitude_kbt=amplitude,sigma=sigma,
                                        flux_ratio=r['common_salt_flux_mol_m2_s']/zero['common_salt_flux_mol_m2_s'],
                                        max_cut_flux_conservation_relative=cut_error,
                                        net_current_relative=r['relative_net_current_error']))
        METRICS['positive_barrier_transport'] = records

    def test_rectangle_continuum_resistance_grid_convergence(self):
        records=[]
        for n in (4,8,16,32):
            g=Geometry(cells_per_fiber=n,fiber_present=False)
            r,y,e,m,h=solve(30.,1.,0.,geometry=g,dp=1.6e-9,dm=1.6e-9,
                            ion_energy_model=CommonIonEnergy('rectangle',3.))
            # Independent 1D exact steady continuum result: electrochemical
            # activity falls across series resistance integral exp(w(x)) dx.
            length_out=2*g.external_fiber_sides_each*g.fiber_side_nm*1e-9
            length_mem=g.fiber_side_nm*1e-9
            exact=1.6e-9*29/(length_out+length_mem*math.exp(3.))
            error=abs(r['common_salt_flux_mol_m2_s']/exact-1)
            records.append(dict(cells_per_fiber=n,flux=r['common_salt_flux_mol_m2_s'],
                                independent_continuum_flux=exact,relative_error=error))
        self.assertTrue(all(a['relative_error']>b['relative_error'] for a,b in zip(records,records[1:])))
        self.assertLess(records[-1]['relative_error'],.015)
        METRICS['rectangle_jump_continuum_grid_convergence'] = records

    def test_no_unvalidated_adsorption_chemistry_mixing(self):
        with self.assertRaisesRegex(ValueError,'fixed sigma only'):
            solve(surface_model={'site_density_nm2':.05},ion_energy_model=CommonIonEnergy(amplitude_kbt=1.))


if __name__=='__main__':
    started=time.perf_counter()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(EnergyTests))
    report=dict(status='PASS' if result.wasSuccessful() else 'FAIL',tests_run=result.testsRun,
                failures=len(result.failures),errors=len(result.errors),runtime_seconds=time.perf_counter()-started,
                metrics=METRICS,scope='prescribed common Na/Cl entry energies; fixed sigma; finite representative 2D water/fiber geometry')
    (ROOT/'energy_test_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    raise SystemExit(0 if result.wasSuccessful() else 1)
