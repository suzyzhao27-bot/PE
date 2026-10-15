"""Small finite 2D staggered five-fiber stack; representative, not reconstruction."""
import csv
import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from solver import Geometry,solve,linear_conductance
from run import CHEM,BRIDGE_SIGMA

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"


def main():
    started=time.perf_counter();rows=[];grid=[]
    for n in (8,16,24):
        g=Geometry(cells_per_fiber=n,fiber_rows=5)
        r0=solve(30,1,0,geometry=g)[0]
        rb,yb,eb,mb,hb=solve(30,1,BRIDGE_SIGMA,geometry=g)
        gb,eq=linear_conductance(1,BRIDGE_SIGMA,geometry=g)
        grid.append(dict(cells_per_fiber=n,grid_nm=rb["grid_nm"],neutral_J=r0["common_salt_flux_mol_m2_s"],
                         bridge_J=rb["common_salt_flux_mol_m2_s"],bridge_neutral_J_ratio=rb["common_salt_flux_mol_m2_s"]/r0["common_salt_flux_mol_m2_s"],
                         bridge_delta_voltage_mv=rb["right_minus_left_theoretical_aqueous_voltage_mv"],
                         bridge_G_s_m2=gb["projected_area_ionic_conductance_s_m2"],bridge_t_plus=gb["linear_current_cation_fraction_t_plus"]))
        if n==24:
            for label,sigma,chem in (("five_rows_neutral",0.,None),("five_rows_HO_reference_frozen",BRIDGE_SIGMA,None),
                                     ("five_rows_hypothetical_NaCl_adsorption",0.,CHEM)):
                r,y,e,m,h=solve(30,1,sigma,geometry=g,surface_model=chem)
                conduct,eq=linear_conductance(1,sigma,geometry=g,surface_model=chem)
                rows.append(dict(scenario=label,transport=r,conductance=conduct,
                                 same_domain_neutral_J_ratio=r["common_salt_flux_mol_m2_s"]/r0["common_salt_flux_mol_m2_s"]))
            np.savez_compressed(OUT/"five_rows_bridge_30_1_fields.npz",mask=mb.mask,ids=mb.ids,
                                u_dimensionless=eb["u"],c_plus_mol_m3=eb["cp"],c_minus_mol_m3=eb["cm"],grid_nm=mb.h)
    with (OUT/"multifiber_grid_study.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(grid[0]));w.writeheader();w.writerows(grid)
    a,b=grid[-2:]
    qa=dict(grid16_24_ratio_difference=abs(a["bridge_neutral_J_ratio"]-b["bridge_neutral_J_ratio"]),
            grid16_24_voltage_difference_mv=abs(a["bridge_delta_voltage_mv"]-b["bridge_delta_voltage_mv"]),
            grid16_24_G_relative_difference=abs(a["bridge_G_s_m2"]-b["bridge_G_s_m2"])/b["bridge_G_s_m2"],
            grid16_24_tplus_difference=abs(a["bridge_t_plus"]-b["bridge_t_plus"]))
    result=dict(status="actual finite 2D five-row representative fibrous geometry; not measured/reconstructed PE network",
                runtime_seconds=time.perf_counter()-started,geometry=asdict(Geometry(cells_per_fiber=24,fiber_rows=5)),
                representative_network_span_nm=239.4,observed_approx_film_thickness_nm=230,
                comparison_is_scale_anchor_not_geometry_identity=True,grid_study=grid,verification=qa,cases=rows,
                code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"solver.py",ROOT/"run_multifiber.py")})
    (OUT/"multifiber_results.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(dict(runtime_seconds=result["runtime_seconds"],verification=qa,
                         cases=[dict(scenario=r["scenario"],J_ratio=r["same_domain_neutral_J_ratio"],
                                     V=r["transport"]["right_minus_left_theoretical_aqueous_voltage_mv"],
                                     G=r["conductance"]["projected_area_ionic_conductance_s_m2"],
                                     t_plus=r["conductance"]["linear_current_cation_fraction_t_plus"]) for r in rows]),indent=2))


if __name__=="__main__":
    main()
