"""Targeted refinement of the primary five-row chemical Robin scenarios."""
import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from solver import Geometry,solve,linear_conductance,RG,F
from run import CHEM

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"


def main():
    started=time.perf_counter();rows=[];cases=[]
    for n in (16,24,32):
        g=Geometry(cells_per_fiber=n,fiber_rows=5)
        r0=solve(30,1,0,geometry=g)[0]
        r,y,e,m,h=solve(30,1,geometry=g,surface_model=CHEM)
        conduct,eq=linear_conductance(1,geometry=g,surface_model=CHEM)
        rows.append(dict(cells_per_fiber=n,grid_nm=r["grid_nm"],J_neutral=r0["common_salt_flux_mol_m2_s"],
                         J_chemical=r["common_salt_flux_mol_m2_s"],chemical_neutral_J_ratio=r["common_salt_flux_mol_m2_s"]/r0["common_salt_flux_mol_m2_s"],
                         delta_voltage_mv=r["right_minus_left_theoretical_aqueous_voltage_mv"],
                         G_s_m2=conduct["projected_area_ionic_conductance_s_m2"],t_plus=conduct["linear_current_cation_fraction_t_plus"],
                         static_sigma_mean_c_m2=conduct["wall_sigma_mean_equilibrium_c_m2"],
                         open_sigma_mean_c_m2=r["wall_sigma_mean_c_m2"],
                         residual=r["max_scaled_equation_residual"],net_current_relative=r["relative_net_current_error"],
                         plus_conservation=r["plus_conservation_relative"],minus_conservation=r["minus_conservation_relative"],
                         Gauss_balance=r["poisson_gauss_balance_scaled"]))
        cases.append(dict(grid=n,neutral_transport=r0,chemical_transport=r,chemical_conductance=conduct))
        if n==32:
            np.savez_compressed(OUT/"five_rows_chemical_30_1_fields.npz",mask=m.mask,ids=m.ids,
                                potential_mv=e["u"]*RG*298.15/F*1000,c_plus_mol_m3=e["cp"],c_minus_mol_m3=e["cm"],
                                sigma_c_m2=e["sigma_local"],grid_nm=m.h)
            es=eq[2];ms=eq[3]
            np.savez_compressed(OUT/"five_rows_chemical_static_1mM_fields.npz",mask=ms.mask,ids=ms.ids,
                                potential_mv=es["u"]*RG*298.15/F*1000,c_plus_mol_m3=es["cp"],c_minus_mol_m3=es["cm"],
                                sigma_c_m2=es["sigma_local"],grid_nm=ms.h)
    a,b=rows[-2:]
    qa=dict(ratio24_32_difference=abs(a["chemical_neutral_J_ratio"]-b["chemical_neutral_J_ratio"]),
            voltage24_32_difference_mv=abs(a["delta_voltage_mv"]-b["delta_voltage_mv"]),
            G24_32_relative_difference=abs(a["G_s_m2"]-b["G_s_m2"])/b["G_s_m2"],
            tplus24_32_difference=abs(a["t_plus"]-b["t_plus"]),
            static_sigma24_32_relative_difference=abs(a["static_sigma_mean_c_m2"]-b["static_sigma_mean_c_m2"])/abs(b["static_sigma_mean_c_m2"]))
    report=dict(status="primary five-row chemical scenarios computed and refined; assumed chemistry/geometry, not PE validation",
                runtime_seconds=time.perf_counter()-started,geometry=asdict(Geometry(cells_per_fiber=32,fiber_rows=5)),
                chemical_model=CHEM,grid_study=rows,verification=qa,cases=cases,
                code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"solver.py",ROOT/"refine_multifiber.py")})
    (OUT/"five_rows_chemical_refinement.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(dict(runtime_seconds=report["runtime_seconds"],finest=rows[-1],verification=qa),indent=2))


if __name__=="__main__":
    main()
