import csv
import hashlib
import json
import platform
import time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from solver import Geometry,solve,linear_conductance,F,RG,SCALE_M

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"
CHEM={"site_density_nm2":.05,"K_Na_per_M":0.,"K_Cl_per_M":10000.}
BRIDGE_SIGMA=-.00674949


def write_csv(name,rows):
    with (OUT/name).open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def selected(row):
    fields=["c_left_mol_m3","c_right_mol_m3","surface_boundary","sigma_c_m2",
            "wall_sigma_mean_c_m2","wall_sigma_min_c_m2","wall_sigma_max_c_m2",
            "plus_flux_right_mol_m2_s","minus_flux_right_mol_m2_s","common_salt_flux_mol_m2_s",
            "right_minus_left_theoretical_aqueous_voltage_mv","relative_net_current_error",
            "plus_conservation_relative","minus_conservation_relative","poisson_gauss_balance_scaled",
            "max_scaled_equation_residual","grid_nm","water_cells","external_water_each_nm"]
    return {key:row[key] for key in fields}


def main():
    OUT.mkdir(exist_ok=True);started=time.perf_counter()
    g=Geometry(cells_per_fiber=32)
    transport=[];all_transport=[]
    for c in (30.,270.):
        neutral=solve(c,1,0,geometry=g)[0]
        for label,sigma,chem in (("fixed_neutral",0.,None),("fixed_negative",-.01,None),
                                 ("fixed_positive",.01,None),("HO_reference_sigma_frozen",BRIDGE_SIGMA,None),
                                 ("hypothetical_NaCl_adsorption_Robin",0.,CHEM)):
            r,y,e,m,h=solve(c,1,sigma,geometry=g,surface_model=chem)
            r["scenario"]=label;r["same_domain_neutral_flux_ratio"]=r["common_salt_flux_mol_m2_s"]/neutral["common_salt_flux_mol_m2_s"]
            row=dict(scenario=label,**selected(r),same_domain_neutral_flux_ratio=r["same_domain_neutral_flux_ratio"])
            transport.append(row);all_transport.append(r)
            if c==30 and label in ("fixed_neutral","fixed_negative","hypothetical_NaCl_adsorption_Robin"):
                np.savez_compressed(OUT/(label+"_30_1_fields.npz"),mask=m.mask,ids=m.ids,
                                    potential_dimensionless=e["u"],c_plus_mol_m3=e["cp"],c_minus_mol_m3=e["cm"],
                                    local_sigma_c_m2=e["sigma_local"],grid_nm=m.h)
    write_csv("transport_cases.csv",transport)
    conductance=[]
    for c in (1.,10.,100.):
        for label,sigma,chem in (("fixed_neutral",0.,None),("fixed_negative",-.01,None),
                                 ("HO_reference_sigma_frozen",BRIDGE_SIGMA,None),
                                 ("hypothetical_NaCl_adsorption_Robin",0.,CHEM)):
            r,eq=linear_conductance(c,sigma,geometry=g,surface_model=chem)
            r["scenario"]=label;conductance.append(r)
    write_csv("linear_G_tplus_static_charge.csv",conductance)
    grid=[]
    for n in (8,16,24,32):
        gg=Geometry(cells_per_fiber=n)
        r0=solve(30,1,0,geometry=gg)[0];r=solve(30,1,BRIDGE_SIGMA,geometry=gg)[0]
        a,eq=linear_conductance(1,BRIDGE_SIGMA,geometry=gg)
        b,eqb=linear_conductance(1,0,geometry=gg,surface_model=CHEM)
        grid.append(dict(cells_per_fiber=n,grid_nm=r["grid_nm"],J_neutral=r0["common_salt_flux_mol_m2_s"],
                         J_charged=r["common_salt_flux_mol_m2_s"],charged_neutral_ratio=r["common_salt_flux_mol_m2_s"]/r0["common_salt_flux_mol_m2_s"],
                         delta_voltage_mv=r["right_minus_left_theoretical_aqueous_voltage_mv"],
                         G=a["projected_area_ionic_conductance_s_m2"],t_plus=a["linear_current_cation_fraction_t_plus"],
                         chemical_sigma_at_1mM=b["wall_sigma_mean_equilibrium_c_m2"],chemical_t_plus=b["linear_current_cation_fraction_t_plus"]))
    write_csv("grid_study.csv",grid)
    domains=[]
    kp=F*F*1.33e-9/(RG*298.15);km=F*F*2.03e-9/(RG*298.15)
    for extent in (2,4,6):
        gg=Geometry(cells_per_fiber=16,external_fiber_sides_each=extent)
        r0=solve(30,1,0,geometry=gg)[0];r=solve(30,1,BRIDGE_SIGMA,geometry=gg)[0]
        a,eq=linear_conductance(1,BRIDGE_SIGMA,geometry=gg)
        total_g=a["projected_area_ionic_conductance_s_m2"];tp=a["linear_current_cation_fraction_t_plus"]
        wallcharge=BRIDGE_SIGMA*4*gg.fiber_side_nm*SCALE_M
        screened=-a["equilibrium_mobile_charge_mol_per_m_depth"]*F/wallcharge
        domains.append(dict(external_water_each_nm=extent*gg.fiber_side_nm,total_height_nm=r["total_aqueous_domain_height_nm"],
                            J_neutral=r0["common_salt_flux_mol_m2_s"],J_charged=r["common_salt_flux_mol_m2_s"],
                            charged_neutral_ratio=r["common_salt_flux_mol_m2_s"]/r0["common_salt_flux_mol_m2_s"],
                            delta_voltage_mv=r["right_minus_left_theoretical_aqueous_voltage_mv"],G=total_g,t_plus=tp,
                            static_screened_wall_charge_fraction=screened,
                            plus_excess_area_resistance_ohm_m2=1/(total_g*tp)-r["total_aqueous_domain_height_nm"]*SCALE_M/kp,
                            minus_excess_area_resistance_ohm_m2=1/(total_g*(1-tp))-r["total_aqueous_domain_height_nm"]*SCALE_M/km,
                            note="Extending axial water adds physical series resistance; total J/G/V/t+ need not be domain invariant. Local excess ionic resistance and static field screening are the domain checks."))
    write_csv("domain_study.csv",domains)
    gridlast,gridprev=grid[-1],grid[-2]
    domainlast,domainprev=domains[-1],domains[-2]
    qa=dict(grid_24_32_ratio_absolute_difference=abs(gridlast["charged_neutral_ratio"]-gridprev["charged_neutral_ratio"]),
            grid_24_32_voltage_difference_mv=abs(gridlast["delta_voltage_mv"]-gridprev["delta_voltage_mv"]),
            grid_24_32_G_relative_difference=abs(gridlast["G"]-gridprev["G"])/gridlast["G"],
            grid_24_32_tplus_absolute_difference=abs(gridlast["t_plus"]-gridprev["t_plus"]),
            domain106_159nm_plus_excess_R_relative_difference=abs(domainlast["plus_excess_area_resistance_ohm_m2"]-domainprev["plus_excess_area_resistance_ohm_m2"])/abs(domainlast["plus_excess_area_resistance_ohm_m2"]),
            domain106_159nm_minus_excess_R_relative_difference=abs(domainlast["minus_excess_area_resistance_ohm_m2"]-domainprev["minus_excess_area_resistance_ohm_m2"])/abs(domainlast["minus_excess_area_resistance_ohm_m2"]),
            static_screening_fraction_in_largest_domain=domainlast["static_screened_wall_charge_fraction"],
            all_transport_max_scaled_residual=max(r["max_scaled_equation_residual"] for r in all_transport),
            all_transport_max_relative_net_current=max(r["relative_net_current_error"] for r in all_transport),
            all_transport_max_relative_ion_flux_conservation=max(max(r["plus_conservation_relative"],r["minus_conservation_relative"]) for r in all_transport))
    (OUT/"verification_summary.json").write_text(json.dumps(qa,indent=2)+"\n")
    manifest=dict(date="2026-10-14",status="actual finite aqueous 2D PNP calculations; hypothetical representative local geometry/chemical parameters",
                  runtime_seconds=time.perf_counter()-started,python=platform.python_version(),geometry=asdict(g),
                  coordinate="positive axial direction from left high-salt to right low-salt boundary; deltaV=V_right-V_left",
                  parameters=dict(T_K=298.15,eps_r=78.3,D_Na_m2_s=1.33e-9,D_Cl_m2_s=2.03e-9,chemical_model=CHEM),
                  bridge_sigma_note="-.00674949 C/m2 is an assumed external H/OH reference surface charge supplied by parent calculation; frozen here, not concentration-dependent H/OH charge regulation",
                  measured_sample_anchor="Polymer2023 draw604 SEM fibre diameter26.6nm and thickness~230nm; square fiber and53.2nm wet gap are idealized, not reconstructed from dry52.4nm pore statistic",
                  scope=dict(no_support=True,finite_entrances_and_external_water=True,local_velocity_everywhere_zero=True,
                             full230nm_membrane_reconstructed=False,real_sigma_or_site_constants_established=False,DI_calculated=False,
                             actual_electrodes_included=False,H_OH_or_pH_chemistry_included=False,three_dimensional_network=False),
                  code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"solver.py",ROOT/"tests.py",ROOT/"run.py")},
                  verification=qa,transport=all_transport,conductance=conductance)
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(dict(runtime_seconds=manifest["runtime_seconds"],transport_cases=len(transport),G_cases=len(conductance),verification=qa),ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
