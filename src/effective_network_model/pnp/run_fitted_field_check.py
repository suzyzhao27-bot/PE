"""Check root-supplied energy coefficients in a different representative 2D geometry.

These are not an absolute fit of the 2D network to experimental curves.
"""
import csv
import json
import time
from pathlib import Path
from solver import CommonIonEnergy, Geometry, solve

ROOT=Path(__file__).resolve().parent


def main():
    start=time.perf_counter();rows=[];details=[]
    for n in (8,16,24):
        g=Geometry(cells_per_fiber=n,fiber_rows=5)
        reference=solve(30.,1.,0.,geometry=g,tol=1e-13)[0]
        amplitudes=(0.,1.,2.,7.03979054) if n==8 else (7.03979054,)
        for amplitude in amplitudes:
            model=CommonIonEnergy('two_sided_exponential',amplitude,893.953885)
            try:
                r,y,e,m,h=solve(30.,1.,0.,geometry=g,tol=1e-13,ion_energy_model=model,max_iter=100)
                row=dict(status='converged',cells_per_fiber=n,ws_kbt=amplitude,ell_nm=model.decay_nm,
                         sigma_given_c_m2=0.,energy_min_water_kbt=r['common_ion_energy_min_kbt'],
                         energy_max_water_kbt=r['common_ion_energy_max_kbt'],
                         flux_mol_m2_s=r['common_salt_flux_mol_m2_s'],
                         flux_ratio_to_same_sigma_zero_energy=r['common_salt_flux_mol_m2_s']/reference['common_salt_flux_mol_m2_s'],
                         voltage_mv=r['right_minus_left_theoretical_aqueous_voltage_mv'],
                         equation_residual=r['max_scaled_equation_residual'],
                         open_current_relative=r['relative_net_current_error'],
                         plus_conservation=r['plus_conservation_relative'],minus_conservation=r['minus_conservation_relative'],
                         gauss_scaled=r['poisson_gauss_balance_scaled'],iterations=r['newton_iterations'],error='')
                details.append(dict(model='root-coefficient transfer spatial check, not geometry-matched fit',result=r))
            except (RuntimeError,ValueError) as error:
                row=dict(status='failed',cells_per_fiber=n,ws_kbt=amplitude,ell_nm=model.decay_nm,
                         sigma_given_c_m2=0.,energy_min_water_kbt=None,energy_max_water_kbt=None,
                         flux_mol_m2_s=None,flux_ratio_to_same_sigma_zero_energy=None,voltage_mv=None,
                         equation_residual=None,open_current_relative=None,plus_conservation=None,
                         minus_conservation=None,gauss_scaled=None,iterations=None,error=str(error))
            rows.append(row)
            print(json.dumps(row),flush=True)
    out=ROOT/'outputs';out.mkdir(exist_ok=True)
    with (out/'root_fitted_energy_spatial_check.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(runtime_seconds=time.perf_counter()-start,requested_root_coefficients=dict(ws_kbt=7.03979054,ell_nm=893.953885),
                 scope='neutral fixed-sigma finite representative five-row square-fiber network; 30/1mM; coefficients transferred for spatial check only',
                 fit_geometry_is_different=True,rows=rows,details=details)
    (out/'root_fitted_energy_spatial_check.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
