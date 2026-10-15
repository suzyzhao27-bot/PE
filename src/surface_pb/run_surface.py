from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.optimize import brentq
from fibre_pb import solve,Adsorption,AcidSites,E,R,F
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs';OUT.mkdir(exist_ok=True)
def clean(v):
 if isinstance(v,np.ndarray):return v.tolist()
 if isinstance(v,np.generic):return v.item()
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
 if isinstance(v,list):return [clean(x) for x in v]
 return v
def slim(r):return {k:clean(v) for k,v in r.items() if k not in ('radius_nm','u_profile','concentration_profiles_mol_m3','species')}
def csvout(name,rows):
 with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
rows=[];profiles=[]
for salt in (0.1,1,10,100,270):
 for sigma in (0,-.001,-.005,-.01,.005):
  r=solve(salt,7,sigma=sigma)
  rows.append(dict(model='fixed_sigma',salt_mM=salt,pH=7,assumed_Gamma_nm2='',assumed_dG_H_kj='',assumed_dG_OH_kj='',sigma_C_m2=r['surface_sigma_C_m2'],surface_mV=r['surface_psi_mV'],slip0p5_mV=r['slip_0p5nm_psi_mV'],edge_mV=r['cell_edge_psi_mV'],K_Na=r['mobile_partition']['Na'],K_Cl=r['mobile_partition']['Cl'],t_cation=r['linear_axial_current_fraction_cations'],t_Na_pair_proxy=r['Na_fraction_NaCl_only_proxy'],gauss_error=r['relative_Gauss_error'],residual=r['residual_scaled']))
for salt in (1,10,100):
 for pH in (3,5,7,9,10):
  for dgoh in (-20,-30,-40,-50,-60):
   a=Adsorption(.05,{'H':-20.,'OH':dgoh});r=solve(salt,pH,adsorption=a)
   rows.append(dict(model='regulated_H_OH_hypothesis',salt_mM=salt,pH=pH,assumed_Gamma_nm2=.05,assumed_dG_H_kj=-20,assumed_dG_OH_kj=dgoh,sigma_C_m2=r['surface_sigma_C_m2'],surface_mV=r['surface_psi_mV'],slip0p5_mV=r['slip_0p5nm_psi_mV'],edge_mV=r['cell_edge_psi_mV'],K_Na=r['mobile_partition']['Na'],K_Cl=r['mobile_partition']['Cl'],t_cation=r['linear_axial_current_fraction_cations'],t_Na_pair_proxy=r['Na_fraction_NaCl_only_proxy'],gauss_error=r['relative_Gauss_error'],residual=r['residual_scaled']))
   if salt==1 and pH in (3,7,10) and dgoh==-50:
    for j,radius in enumerate(r['radius_nm']):
     profiles.append(dict(pH=pH,salt_mM=salt,r_nm=radius,u=r['u_profile'][j],psi_mV=r['u_profile'][j]*R*298.15/F*1000,**{n+'_mol_m3':r['concentration_profiles_mol_m3'][k,j] for k,n in enumerate(r['species'])}))
for row in rows:
 base=solve(row['salt_mM'],row['pH'],sigma=0)
 row['same_water_uncharged_t_cation']=base['linear_axial_current_fraction_cations']
 row['delta_t_cation_from_same_water']=row['t_cation']-row['same_water_uncharged_t_cation']
 row['relative_Na_Cl_partition_ratio']=row['K_Na']/row['K_Cl']
csvout('surface_scenarios.csv',rows);csvout('surface_radial_profiles.csv',profiles)
thresholds=[]
# This inverse calculation targets a MODEL surface potential, never measured zeta.
for salt in (0.1,1,10,100):
 target=-70.
 sig=brentq(lambda sigma:solve(salt,7,sigma=sigma)['surface_psi_mV']-target,-.15,-1e-7,xtol=1e-13)
 gammin=-sig/(E*1e18)
 threshold=dict(salt_mM=salt,target_assumed_surface_mV=target,required_fixed_sigma_C_m2=sig,minimum_full_occupancy_sites_nm2=gammin,assumed_sites_nm2=.05)
 if gammin>=.05:
  threshold.update(OH_binding_dG_kj_for_Gamma0p05=None,status='impossible at declared site density even at full adsorption')
 else:
  dg=brentq(lambda dg:solve(salt,7,adsorption=Adsorption(.05,{'OH':dg}))['surface_psi_mV']-target,-100,0,xtol=1e-7)
  threshold.update(OH_binding_dG_kj_for_Gamma0p05=dg,status='hypothetical affinity threshold, not measured PE energy')
 thresholds.append(threshold)
acid_rows=[]
for salt in (1,10,100):
 for pH in (3,7,10):
  for pk in (4,7,10):
   acid=solve(salt,pH,adsorption=AcidSites(.1,pk));acid_rows.append({'salt_mM':salt,'pH':pH,'assumed_pKa':pk,'Gamma_nm2':.1,'sigma_C_m2':acid['surface_sigma_C_m2'],'psi_surface_mV':acid['surface_psi_mV'],'deprotonated_fraction':acid['surface_adsorption_occupancies']['deprotonated'],'t_cation':acid['linear_axial_current_fraction_cations'],'gauss_error':acid['relative_Gauss_error']})
csvout('acid_site_hypotheses.csv',acid_rows)
(OUT/'surface_thresholds.json').write_text(json.dumps(clean(thresholds),ensure_ascii=False,indent=2))
# Selected origin demonstration, full machine-readable profiles and all controls.
selected=[]
for pH in (3,7,10):selected.append(slim(solve(1,pH,adsorption=Adsorption(.05,{'H':-20,'OH':-50}))))
manifest={'scope':'exterior solid fibre PB plus self-consistent bound-ion adsorption; no imposed sigma in regulated runs','geometry_source':'Polymer2023 draw604 SEM mean d26.6nm; cell b39.5nm/fully connected water are assumptions','hypotheses':{'Gamma_nm2':.05,'dG_H_kJ_mol':-20,'dG_OH_kJ_mol':-50,'free_energies_not_PE_measurements':True,'slip_plane_nm_assumed':.5,'ideal_point_ions_no_flow':True},'selected_regulated':selected,'scenario_rows':len(rows),'profiles_rows':len(profiles),'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}}
(OUT/'surface_manifest.json').write_text(json.dumps(clean(manifest),ensure_ascii=False,indent=2))
print(json.dumps({'selected':[{k:x[k] for k in ('pH','surface_sigma_C_m2','surface_psi_mV','linear_axial_current_fraction_cations','relative_Gauss_error')} for x in selected],'thresholds':thresholds,'rows':len(rows)},ensure_ascii=False,indent=2))
