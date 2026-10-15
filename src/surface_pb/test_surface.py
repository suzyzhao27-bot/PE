import unittest,json,math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp
from fibre_pb import solve,Adsorption,AcidSites,reservoir,F,R,EPS0,E,Z
METRICS={}
class Physics(unittest.TestCase):
 def test_neutral_and_chemical_pH_baseline(self):
  for pH in (3,7,10):
   v=solve(1,pH,sigma=0);self.assertLess(abs(v['surface_psi_mV']),1e-9)
   self.assertAlmostEqual(v['linear_axial_current_fraction_cations'],v['same_reservoir_uncharged_cation_fraction'],places=11)
   self.assertLess(abs(sum(Z[k]*x for k,x in reservoir(1,pH).items())),1e-10)
  METRICS['pH3_neutral_current_cation_fraction']=solve(1,3,sigma=0)['linear_axial_current_fraction_cations']
 def test_cell_Gauss_and_charge_sign(self):
  for sig in (-.01,-.001,.001,.01):
   v=solve(1,7,sigma=sig);self.assertLess(v['relative_Gauss_error'],1e-8)
   self.assertGreater(v['surface_psi_mV']*sig,0)
   self.assertLess(v['surface_bound_charge_per_2pi_length_C_m']*v['mobile_charge_per_2pi_length_C_m'],0)
 def test_fixed_charge_symmetry(self):
  a=solve(1,sigma=-.005);b=solve(1,sigma=.005)
  self.assertAlmostEqual(a['surface_psi_mV'],-b['surface_psi_mV'],places=8)
  self.assertAlmostEqual(a['mobile_partition']['Na'],b['mobile_partition']['Cl'],places=9)
 def test_regulated_zero_affinity_symmetry_and_capacity(self):
  sym=solve(1,7,adsorption=Adsorption(.05,{'Na':0.,'Cl':0.}))
  self.assertLess(abs(sym['surface_sigma_C_m2']),1e-12)
  for pH in (3,7,10):
   v=solve(1,pH,adsorption=Adsorption(.05,{'H':-20,'OH':-50}))
   self.assertLessEqual(abs(v['surface_sigma_C_m2']),E*.05*1e18*(1+1e-12))
   self.assertLess(v['relative_Gauss_error'],1e-8)
 def test_grid_refinement(self):
  vals=[solve(1,7,n=n,adsorption=Adsorption(.05,{'H':-20,'OH':-50})) for n in (128,256,512,1024)]
  psi=[v['surface_psi_mV'] for v in vals];delta=np.abs(np.diff(psi))
  self.assertLess(delta[-1],delta[-2]*.35);self.assertLess(delta[-1],.01)
  METRICS['grid']={'n':[128,256,512,1024],'surface_mV':psi,'last_mV_difference':delta[-1]}
 def test_independent_collocation_solver(self):
  bulk=reservoir(1,7);I=bulk['Na']+bulk['H'];chi=F*F/(EPS0*78.3*R*298.15)*1e-18
  a,b=13.3,39.5;gamma=.05*1e18;dgh,dgo=-20.,-50.
  def sc(u):
   bh=math.exp(-dgh*1000/(R*298.15))*bulk['H']/1000*math.exp(-u)
   bo=math.exp(-dgo*1000/(R*298.15))*bulk['OH']/1000*math.exp(u)
   return E*gamma*(bh-bo)/(1+bh+bo)
  def eq(r,y):return np.vstack((y[1],-y[1]/r+2*chi*I*np.sinh(y[0])))
  def bc(ya,yb):return np.array([ya[1]+F/(EPS0*78.3*R*298.15)*1e-9*sc(float(ya[0])),yb[1]])
  fv=solve(1,7,n=2048,adsorption=Adsorption(.05,{'H':-20,'OH':-50}))
  mesh=np.linspace(a,b,180);u=np.interp(mesh,fv['radius_nm'],fv['u_profile']);y=np.vstack((u,np.gradient(u,mesh)))
  sol=solve_bvp(eq,bc,mesh,y,tol=1e-9,max_nodes=10000)
  self.assertTrue(sol.success,sol.message)
  bp=sol.sol(a)[0]*R*298.15/F*1000;bs=sc(float(sol.sol(a)[0]))
  self.assertLess(abs(bp-fv['surface_psi_mV']),.002)
  self.assertLess(abs(bs-fv['surface_sigma_C_m2']),1e-6)
  METRICS['independent_BVP']={'surface_mV':bp,'FV2048_surface_mV':fv['surface_psi_mV'],'absolute_difference_mV':abs(bp-fv['surface_psi_mV']),'sigma_BVP':bs,'sigma_FV':fv['surface_sigma_C_m2'],'bvp_max_residual':float(sol.rms_residuals.max())}
 def test_planar_Grahame_limit(self):
  sigma=-.005;c=reservoir(1,7)['Na']+reservoir(1,7)['H'];eps=EPS0*78.3
  kappa=math.sqrt(2*F*F*c/(eps*R*298.15))
  exact=2*math.asinh(sigma*F/(2*eps*R*298.15*kappa))*R*298.15/F*1000
  v=solve(1,7,radius_nm=100000,cell_nm=100150,n=2048,sigma=sigma)
  self.assertLess(abs(v['surface_psi_mV']-exact),.02)
  METRICS['planar_limit']={'model_mV':v['surface_psi_mV'],'Grahame_mV':exact,'difference_mV':abs(v['surface_psi_mV']-exact)}
 def test_outer_bulk_flux_and_domain(self):
  v=solve(1,7,radius_nm=13.3,cell_nm=130,n=2048,sigma=-.005,outer='bulk')
  self.assertLess(v['relative_Gauss_error'],1e-8);self.assertGreater(abs(v['outer_electric_flux_C_m']),0)
  a=solve(1,7,radius_nm=13.3,cell_nm=130,n=1024,sigma=-.005,outer='bulk')
  b=solve(1,7,radius_nm=13.3,cell_nm=180,n=1536,sigma=-.005,outer='bulk')
  self.assertLess(abs(a['surface_psi_mV']-b['surface_psi_mV']),.01)
  METRICS['bulk_outer_domain']={'r130_mV':a['surface_psi_mV'],'r180_mV':b['surface_psi_mV'],'difference_mV':abs(a['surface_psi_mV']-b['surface_psi_mV'])}
 def test_acid_site_direction_and_pKa(self):
  v=[solve(1,pH,adsorption=AcidSites(.1,7)) for pH in (3,7,10)]
  self.assertTrue(abs(v[0]['surface_sigma_C_m2'])<abs(v[1]['surface_sigma_C_m2'])<abs(v[2]['surface_sigma_C_m2']))
  weak=solve(1,7,adsorption=AcidSites(.1,10))
  METRICS['declared_pKa10_pH7']={'sigma':weak['surface_sigma_C_m2'],'surface_mV':weak['surface_psi_mV']}
if __name__=='__main__':
 result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Physics))
 def convert(v):
  if isinstance(v,np.generic):return v.item()
  if isinstance(v,dict):return {k:convert(x) for k,x in v.items()}
  return v
 qa_path=Path(__file__).resolve().parent.parent/'qa/surface_tests.json'
 qa_path.parent.mkdir(parents=True,exist_ok=True)
 qa_path.write_text(json.dumps(convert({'status':'PASS' if result.wasSuccessful() else 'FAIL','tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'metrics':METRICS}),ensure_ascii=False,indent=2))
 raise SystemExit(not result.wasSuccessful())
