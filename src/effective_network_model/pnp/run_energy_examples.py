"""Small representative runs and grid study, no parameter fit to experiments."""
import csv
import hashlib
import json
import platform
import time
from dataclasses import asdict
from pathlib import Path
import numpy as np
import scipy
from solver import Geometry,CommonIonEnergy,solve

ROOT=Path(__file__).resolve().parent


def main():
    started=time.perf_counter();out=ROOT/'outputs';out.mkdir(exist_ok=True)
    rows=[];full=[]
    for n in (8,16,24):
        g=Geometry(cells_per_fiber=n,fiber_rows=5)
        for sigma in (0.,-.005):
            baseline=solve(30.,1.,sigma,geometry=g)[0]
            for label,model in (('zero',None),('rectangle_w3',CommonIonEnergy('rectangle',3.,26.6)),
                                ('two_sided_ws3_ell26p6',CommonIonEnergy('two_sided_exponential',3.,26.6))):
                r,y,e,m,h=solve(30.,1.,sigma,geometry=g,ion_energy_model=model) if model else (baseline,None,None,None,None)
                row=dict(scenario=label,cells_per_fiber=n,sigma_given_c_m2=sigma,
                         network_span_nm=r['representative_network_span_nm'],fiber_square_side_nm=g.fiber_side_nm,
                         nominal_periodic_open_gap_nm=r['nominal_open_gap_nm'],
                         J_mol_m2_s=r['common_salt_flux_mol_m2_s'],
                         J_over_same_sigma_zero_energy=r['common_salt_flux_mol_m2_s']/baseline['common_salt_flux_mol_m2_s'],
                         voltage_mv=r['right_minus_left_theoretical_aqueous_voltage_mv'],
                         max_residual=r['max_scaled_equation_residual'],current_relative=r['relative_net_current_error'],
                         plus_conservation=r['plus_conservation_relative'],minus_conservation=r['minus_conservation_relative'],
                         gauss_scaled=r['poisson_gauss_balance_scaled'])
                rows.append(row);full.append(dict(scenario=label,cells_per_fiber=n,geometry=asdict(g),result=r))
                if n==16 and sigma==0. and model:
                    np.savez_compressed(out/(label+'_example_fields.npz'),mask=m.mask,ids=m.ids,
                                        potential_dimensionless=e['u'],c_plus_mol_m3=e['cp'],c_minus_mol_m3=e['cm'],
                                        common_energy_kbt=e['common_energy_kbt'],grid_nm=m.h)
    with (out/'representative_energy_transport.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    grid=[]
    for sigma in (0.,-.005):
        for label in ('rectangle_w3','two_sided_ws3_ell26p6'):
            selected=[r for r in rows if r['scenario']==label and r['sigma_given_c_m2']==sigma]
            prev,last=selected[-2:]
            grid.append(dict(scenario=label,sigma_given_c_m2=sigma,grid_pair=[16,24],
                             flux_relative_difference=abs(last['J_mol_m2_s']/prev['J_mol_m2_s']-1),
                             ratio_absolute_difference=abs(last['J_over_same_sigma_zero_energy']-prev['J_over_same_sigma_zero_energy']),
                             voltage_difference_mv=abs(last['voltage_mv']-prev['voltage_mv'])))
    report=dict(status='representative calculations completed',runtime_seconds=time.perf_counter()-started,
                python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                parameters=dict(c_left_mM=30,c_right_mM=1,T_K=298.15,eps_r=78.3,Dplus_m2_s=1.33e-9,Dminus_m2_s=2.03e-9),
                energy_definition='w=U/kBT shared by Na and Cl, imposed effective structure/interface entry energy; not dielectric microscopic self-energy',
                geometry_boundary='field applies in network span between first/last square-fiber exterior faces; exterior baths zero; SG uses Delta(z*u+w)',
                geometry_warning='26.6nm square-fiber side and 53.2nm periodic open gap are representative assumptions, not the reduced straight-pore radius a or a 3D reconstruction',
                chemical_model='fixed given sigma only; no old gamma/K adsorption-chemistry mixture',
                grid_study=grid,results=full,
                source_solver_sha256=hashlib.sha256((ROOT/'solver_v3_baseline.py').read_bytes()).hexdigest(),
                code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'solver.py',ROOT/'tests_energy.py',ROOT/'run_energy_examples.py')})
    (out/'representative_energy_results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(runtime_seconds=report['runtime_seconds'],cases=len(rows),grid_study=grid),indent=2))


if __name__=='__main__':main()
