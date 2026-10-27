"""Run the frozen v7 C3 cylinder transport and equilibrium-capacity matrix."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from v7_axisymmetric_pnp import Cylinder, D_CL, D_NA, solve


ROOT = Path(__file__).resolve().parents[2]
SPEC = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json"
OUT = ROOT / "WP2"


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    c3 = json.loads(SPEC.read_text(encoding="utf-8-sig"))["matrices"]["C3"]
    nr,nz = c3["radial_axial_grid_main"]
    transport, capacity = [], []
    for a_nm in c3["a_nm"]:
        cylinder = Cylinder(a_nm*1e-9,a_nm*1e-9*c3["length_over_radius"],nr,nz)
        for c_r in c3["c_R_mol_m3"]:
            for sigma_mc in c3["sigma_mC_m2"]:
                sigma = sigma_mc*1e-3
                for ratio in c3["c_L_over_c_R"]:
                    try:
                        result,_,_ = solve(ratio*c_r,c_r,sigma,cylinder,
                                           c3["T_K"],c3["epsilon_r"])
                        assert result["current_relative"] <= c3["relative_current_tolerance"]
                        assert result["min_c_mol_m3"] > 0
                        assert result["plus_conservation_relative"] <= 1e-6
                        assert result["minus_conservation_relative"] <= 1e-6
                        row = {"a_nm":a_nm,"c_R_mol_m3":c_r,"sigma_mC_m2":sigma_mc,
                               "c_L_over_c_R":ratio, **result, "status":"PASS","error":""}
                    except Exception as exc:
                        row = {"a_nm":a_nm,"c_R_mol_m3":c_r,"sigma_mC_m2":sigma_mc,
                               "c_L_over_c_R":ratio,"status":"FAIL","error":repr(exc)}
                    transport.append(row)
                    if len(transport)%6==0:
                        print(f"C3 transport {len(transport)}/36; failures "
                              f"{sum(r['status']=='FAIL' for r in transport)}", flush=True)
                try:
                    center,_,_ = solve(c_r,c_r,sigma,cylinder,c3["T_K"],c3["epsilon_r"])
                    derivatives = {}
                    for step in c3["equilibrium_capacity_relative_steps"]:
                        plus,_,_ = solve(c_r*(1+step),c_r*(1+step),sigma,cylinder,
                                         c3["T_K"],c3["epsilon_r"])
                        minus,_,_ = solve(c_r*(1-step),c_r*(1-step),sigma,cylinder,
                                          c3["T_K"],c3["epsilon_r"])
                        derivatives[step] = ((plus["mean_salt_concentration_mol_m3"]-
                                              minus["mean_salt_concentration_mol_m3"])/
                                             (2*step*c_r))
                    steps = c3["equilibrium_capacity_relative_steps"]
                    k_inc = derivatives[steps[0]]
                    diff = abs(k_inc-derivatives[steps[1]])/max(abs(derivatives[steps[1]]),1e-300)
                    assert diff <= c3["conductance_capacity_change_tolerance"]
                    assert k_inc > 0
                    row = {"a_nm":a_nm,"c_b_mol_m3":c_r,"sigma_mC_m2":sigma_mc,
                           "nr":nr,"nz":nz,"mean_salt_at_center_mol_m3":
                           center["mean_salt_concentration_mol_m3"],
                           "K_inc_step_1pct":k_inc,"K_inc_step_halfpct":derivatives[steps[1]],
                           "step_relative_difference":diff,"status":"PASS","error":""}
                except Exception as exc:
                    row = {"a_nm":a_nm,"c_b_mol_m3":c_r,"sigma_mC_m2":sigma_mc,
                           "nr":nr,"nz":nz,"status":"FAIL","error":repr(exc)}
                capacity.append(row)
                if len(capacity)%6==0:
                    print(f"C3 capacity {len(capacity)}/18; failures "
                          f"{sum(r['status']=='FAIL' for r in capacity)}", flush=True)
    assert len(transport)==c3["transport_cases"]
    assert len(capacity)==c3["equilibrium_parameter_sets"]
    # Supply stable column sets even when a failed row occurs before any PASS.
    tcols = ["a_nm","c_R_mol_m3","sigma_mC_m2","c_L_over_c_R"]
    ccols = ["a_nm","c_b_mol_m3","sigma_mC_m2","nr","nz"]
    for rows, cols in ((transport,tcols),(capacity,ccols)):
        for row in rows:
            for key in row:
                if key not in cols:
                    cols.append(key)
    for row in transport:
        for key in tcols:
            row.setdefault(key,"")
    for row in capacity:
        for key in ccols:
            row.setdefault(key,"")
    with (OUT/"v7_c3_transport.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=tcols);writer.writeheader();writer.writerows(transport)
    with (OUT/"v7_c3_capacity.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=ccols);writer.writeheader();writer.writerows(capacity)
    failures = sum(r["status"]=="FAIL" for r in transport+capacity)
    print(f"C3 transport={len(transport)}, capacity={len(capacity)}, failures={failures}")
    if failures:
        raise RuntimeError(f"C3 has {failures} failures; see CSV")


if __name__=="__main__":
    main()
