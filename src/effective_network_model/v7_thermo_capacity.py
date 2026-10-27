"""v7 fixed-geometry salt capacity and focused closed-system free-energy checks."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from v7_network_transport import SPEC_PATH, closed_baths, system


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP2"


def free_energy(states: np.ndarray, capacities: np.ndarray) -> np.ndarray:
    """Dimensionless F/(RT C0), with 0 log 0 taken by continuity."""
    positive = np.where(states > 0, states, 1.0)
    terms = np.where(states > 0, states * (np.log(positive) - 1), 0.0)
    return capacities @ terms


def main() -> None:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    config_hash = hashlib.sha256(SPEC_PATH.read_bytes()).hexdigest()
    a_ref = spec["physics"]["a_ref_m"]
    rows = []
    selected = None
    with gzip.open(OUT / "v7_network_tables.jsonl.gz", "rt", encoding="utf-8") as stream:
        for line in stream:
            graph = json.loads(line)
            if (graph["series"], graph["h_nm"], graph["seed"]) == ("H", 200, 7001):
                selected = graph
            xs = sorted({node["x_m"] for node in graph["nodes"]})
            ys = sorted({node["y_m"] for node in graph["nodes"]})
            s = xs[1] - xs[0]
            area = len(xs) * len(ys) * s*s
            h = graph["h_nm"]*1e-9
            water_volume = sum(node["V_m3"] for node in graph["nodes"])
            for beta in (0.0,3.0):
                capacity = sum(node["V_m3"]*np.exp(-beta*a_ref/node["r_m"]) for node in graph["nodes"])
                assert capacity <= water_volume*(1+1e-12)
                assert water_volume <= area*h*(1+1e-12)
                rows.append({"config_sha256":config_hash,"series":graph["series"],"h_nm":graph["h_nm"],
                             "seed":graph["seed"],"beta":beta,"area_m2":area,
                             "water_volume_m3":water_volume,"salt_capacity_m3":capacity,
                             "water_over_Ah":water_volume/(area*h),
                             "salt_capacity_over_Ah":capacity/(area*h),
                             "alpha1000_salt_capacity_over_VR":1000*capacity/area})
    assert len(rows)==192 and selected is not None
    with (OUT/"v7_geometry_capacity_bound.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    results=[]
    for beta in (0.0,3.0):
        capacities,edges,_area=system(selected,beta,1000,spec)
        for state in ("main_zero", "positive"):
            initial=np.array([1.0,0.0]+[0.0]*(len(capacities)-2))
            if state=="positive":
                initial[1]=0.5
                initial[2:]=0.25
            times=np.array([0,0.001,0.01,0.1,1,10,60,120])*60.0
            run=closed_baths(capacities,edges,times,initial=initial,rtol=1e-9,atol=1e-13)
            f=free_energy(run["all_u"],capacities)
            scale=max(abs(f[0]),float(capacities@initial))
            drift=float(np.max(np.diff(f)/scale))
            assert drift <= 1e-10, (beta,state,drift)
            assert run["mass_relative_error"] <= 1e-8
            assert run["minimum_u"] >= -1e-12
            results.append({"beta":beta,"initial":state,"times_min":(times/60).tolist(),
                            "normalized_F":(f/scale).tolist(),"max_upward_F_step_over_scale":drift,
                            "mass_relative_error":run["mass_relative_error"],"minimum_u":run["minimum_u"],
                            "normalization":"max(abs(F_initial),initial_total_salt); F/(RT*C0) in volume units"})
    (OUT/"v7_free_energy_checks.json").write_text(json.dumps({"config_sha256":config_hash,"runs":results},indent=2)+"\n",encoding="utf-8")
    print(f"capacity rows={len(rows)}; max water/Ah={max(r['water_over_Ah'] for r in rows):.6g}; "
          f"max salt/Ah={max(r['salt_capacity_over_Ah'] for r in rows):.6g}")
    print(f"T5 focused free energy runs={len(results)} PASS; max upward normalized step="
          f"{max(r['max_upward_F_step_over_scale'] for r in results):.3g}")


if __name__=="__main__":
    main()
