"""Compare cylinder small-drive conductance and incremental capacity to NS."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP2"


def read(name):
    with (OUT/name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    transport = read("v7_c3_transport.csv")
    capacity = read("v7_c3_capacity.csv")
    assert len(transport)==36 and len(capacity)==18
    assert all(r["status"]=="PASS" for r in transport+capacity)
    small = {(int(r["a_nm"]),int(r["c_R_mol_m3"]),float(r["sigma_mC_m2"])):r
             for r in transport if float(r["c_L_over_c_R"])==1.01}
    neutral = {(int(r["a_nm"]),int(r["c_R_mol_m3"])):
               float(r["F_salt_mol_s"])/(0.01*float(r["c_R_mol_m3"]))
               for r in transport if float(r["sigma_mC_m2"])==0 and float(r["c_L_over_c_R"])==1.01}
    rows = []
    nearest = []
    for cap in capacity:
        a, c, sigma = int(cap["a_nm"]),int(cap["c_b_mol_m3"]),float(cap["sigma_mC_m2"])
        flux = float(small[a,c,sigma]["F_salt_mol_s"])
        g_ratio = flux/(0.01*c*neutral[a,c])
        k_inc = float(cap["K_inc_step_1pct"])
        point = []
        for beta in (0.0,1.0,3.0,6.0):
            k_ns = math.exp(-beta*25/a)
            g_error = abs(k_ns-g_ratio)/max(abs(g_ratio),1e-300)
            k_error = abs(k_ns-k_inc)/max(abs(k_inc),1e-300)
            row = {"a_nm":a,"c_b_mol_m3":c,"sigma_mC_m2":sigma,"beta":beta,
                   "PNP_conductance_over_neutral":g_ratio,"PNP_K_inc":k_inc,
                   "NS_shared_K_and_g_ratio":k_ns,
                   "conductance_relative_error":g_error,"capacity_relative_error":k_error,
                   "joint_max_relative_error":max(g_error,k_error),
                   "within_20pct_both":g_error<=0.2 and k_error<=0.2}
            rows.append(row);point.append(row)
        nearest.append(min(point,key=lambda r:(r["joint_max_relative_error"],r["beta"])))
    assert len(rows)==72 and len(nearest)==18
    worst=max(nearest,key=lambda r:(r["joint_max_relative_error"],-r["a_nm"],-r["c_b_mol_m3"],-r["sigma_mC_m2"]))
    with (OUT/"v7_c3_ns_comparison.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (OUT/"v7_c3_nearest_diagnostic.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(nearest[0]));writer.writeheader();writer.writerows(nearest)
    print(f"C3 fixed-beta rows={len(rows)}; nearest diagnostic within 20% both="
          f"{sum(r['within_20pct_both'] for r in nearest)}/18; worst point "
          f"a={worst['a_nm']} nm c={worst['c_b_mol_m3']} mol/m3 sigma={worst['sigma_mC_m2']} mC/m2 "
          f"joint error={worst['joint_max_relative_error']:.3g}")


if __name__=="__main__":
    main()
