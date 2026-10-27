"""Prescribed 48x480 C3 refinement: four large drives and one worst point."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from v7_axisymmetric_pnp import Cylinder, solve


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"WP2"
SPEC=ROOT/"原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json"


def read(name):
    with (OUT/name).open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    c3=json.loads(SPEC.read_text(encoding="utf-8-sig"))["matrices"]["C3"]
    main_transport=read("v7_c3_transport.csv")
    main_capacity=read("v7_c3_capacity.csv")
    nearest=read("v7_c3_nearest_diagnostic.csv")
    worst=max(nearest,key=lambda r:(float(r["joint_max_relative_error"]),
                                    -int(r["a_nm"]),-int(r["c_b_mol_m3"]),
                                    -float(r["sigma_mC_m2"])))
    cases=[(a,c,-1.0,10.0) for a in (5,50) for c in (1,100)]
    cases.append((int(worst["a_nm"]),int(worst["c_b_mol_m3"]),
                  float(worst["sigma_mC_m2"]),1.01))
    assert len(cases)==5 and len(set(cases))==5
    nr,nz=c3["refinement"]
    rows=[]
    for a,c,sigma_mc,ratio in cases:
        cylinder=Cylinder(a*1e-9,a*1e-9*c3["length_over_radius"],nr,nz)
        old=next(r for r in main_transport if int(r["a_nm"])==a and
                 int(r["c_R_mol_m3"])==c and float(r["sigma_mC_m2"])==sigma_mc and
                 float(r["c_L_over_c_R"])==ratio)
        try:
            result,_,_=solve(c*ratio,c,sigma_mc*1e-3,cylinder,c3["T_K"],c3["epsilon_r"])
            relative=abs(result["F_salt_mol_s"]-float(old["F_salt_mol_s"]))/max(
                abs(result["F_salt_mol_s"]),1e-300)
            status="PASS" if relative<=c3["conductance_capacity_change_tolerance"] and \
                result["current_relative"]<=c3["relative_current_tolerance"] else "NOT_MET"
            row={"kind":"transport","a_nm":a,"c_R_mol_m3":c,
                 "sigma_mC_m2":sigma_mc,"ratio":ratio,"nr":nr,"nz":nz,
                 "coarse_value":float(old["F_salt_mol_s"]),"fine_value":result["F_salt_mol_s"],
                 "relative_change":relative,"current_relative":result["current_relative"],
                 "plus_conservation_relative":result["plus_conservation_relative"],
                 "minus_conservation_relative":result["minus_conservation_relative"],
                 "status":status,"error":""}
        except Exception as exc:
            row={"kind":"transport","a_nm":a,"c_R_mol_m3":c,
                 "sigma_mC_m2":sigma_mc,"ratio":ratio,"nr":nr,"nz":nz,
                 "status":"FAIL","error":repr(exc)}
        rows.append(row)
        print(f"C3 refinement transport {sum(r['kind']=='transport' for r in rows)}/5: {row['status']}",flush=True)
    unique=sorted(set((a,c,sigma_mc) for a,c,sigma_mc,_ in cases))
    for a,c,sigma_mc in unique:
        cylinder=Cylinder(a*1e-9,a*1e-9*c3["length_over_radius"],nr,nz)
        old=next(r for r in main_capacity if int(r["a_nm"])==a and
                 int(r["c_b_mol_m3"])==c and float(r["sigma_mC_m2"])==sigma_mc)
        try:
            step=c3["equilibrium_capacity_relative_steps"][0]
            plus,_,_=solve(c*(1+step),c*(1+step),sigma_mc*1e-3,cylinder,
                            c3["T_K"],c3["epsilon_r"])
            minus,_,_=solve(c*(1-step),c*(1-step),sigma_mc*1e-3,cylinder,
                             c3["T_K"],c3["epsilon_r"])
            k=(plus["mean_salt_concentration_mol_m3"]-
               minus["mean_salt_concentration_mol_m3"])/(2*step*c)
            relative=abs(k-float(old["K_inc_step_1pct"]))/max(abs(k),1e-300)
            row={"kind":"capacity","a_nm":a,"c_R_mol_m3":c,
                 "sigma_mC_m2":sigma_mc,"ratio":1.0,"nr":nr,"nz":nz,
                 "coarse_value":float(old["K_inc_step_1pct"]),"fine_value":k,
                 "relative_change":relative,"current_relative":"",
                 "plus_conservation_relative":"","minus_conservation_relative":"",
                 "status":"PASS" if relative<=c3["conductance_capacity_change_tolerance"] else "NOT_MET",
                 "error":""}
        except Exception as exc:
            row={"kind":"capacity","a_nm":a,"c_R_mol_m3":c,
                 "sigma_mC_m2":sigma_mc,"ratio":1.0,"nr":nr,"nz":nz,
                 "status":"FAIL","error":repr(exc)}
        rows.append(row)
        print(f"C3 refinement capacity {sum(r['kind']=='capacity' for r in rows)}/{len(unique)}: {row['status']}",flush=True)
    assert len(rows)==5+len(unique)
    columns=[]
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    for row in rows:
        for key in columns:
            row.setdefault(key,"")
    with (OUT/"v7_c3_refinement.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    bad=sum(r["status"]!="PASS" for r in rows)
    print(f"C3 48x480 rows={len(rows)}; not passed={bad}")
    if bad:
        raise RuntimeError(f"C3 refinement has {bad} unresolved rows; see CSV")


if __name__=="__main__":
    main()
