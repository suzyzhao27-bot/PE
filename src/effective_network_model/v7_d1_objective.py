"""Frozen v7 D1 B-series objective and deterministic starts; no search here."""

from __future__ import annotations

import csv
import gzip
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from v7_network_transport import SPEC_PATH, closed_baths, system


ROOT=Path(__file__).resolve().parents[2]
TASK=ROOT/"原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0"
OUT=ROOT/"WP4"
H_VALUES=(65,120,200,350,700,2600)
_WORKER_SPEC=None
_WORKER_CURVES=None
_WORKER_GRAPHS=None


def curves():
    roles={}
    with (TASK/"inputs/salt/curve_roles_v7.csv").open(encoding="utf-8-sig",newline="") as stream:
        for row in csv.DictReader(stream):
            roles[row["curve_id"]]=row
    by_h={h:[] for h in H_VALUES}
    with (TASK/"inputs/salt/salt_transport_tidy.csv").open(encoding="utf-8-sig",newline="") as stream:
        for row in csv.DictReader(stream):
            if roles[row["curve_id"]]["role"]!="B":
                continue
            assert row["time_unit"]=="min" and row["conductivity_unit"]=="μS/cm"
            h=int(row["sample_id"].split()[0])
            by_h[h].append((float(row["time"]),float(row["conductivity"]),row["curve_id"]))
    for h,points in by_h.items():
        points.sort()
        assert len(points)==7 and points[0][0]==0
        assert len({t for t,_,_ in points})==7
    assert sum(len(p)-1 for p in by_h.values())==36
    return by_h


def graphs():
    by_h={h:[] for h in H_VALUES}
    with gzip.open(ROOT/"WP2/v7_network_tables.jsonl.gz","rt",encoding="utf-8") as stream:
        for line in stream:
            graph=json.loads(line)
            if graph["series"]=="S":
                by_h[graph["h_nm"]].append(graph)
    for h, group in by_h.items():
        group.sort(key=lambda g:g["seed"])
        assert [g["seed"] for g in group]==list(range(7001,7009))
    return by_h


def starts(model: str):
    assert model in ("N0","NS")
    bounds=[(0.0,5.0),(0.0,6.0)] + ([(0.0,6.0)] if model=="NS" else [])
    mid=[(lo+hi)/2 for lo,hi in bounds]
    rng=np.random.Generator(np.random.PCG64(20261022))
    # Eleven stratified permuted positions in every transformed coordinate.
    jitter=np.column_stack([(rng.permutation(11)+rng.random(11))/11 for _ in bounds])
    mapped=np.array([lo+(hi-lo)*jitter[:,k] for k,(lo,hi) in enumerate(bounds)]).T
    return np.vstack([mid,mapped])


def predict(curves_by_h, graphs_by_h, spec, alpha: float, beta: float, q: float):
    predictions={}
    diagnostics={}
    for h in H_VALUES:
        times=np.array([point[0] for point in curves_by_h[h]])*60.0
        right=[]
        mass_errors=[]
        minima=[]
        for graph in graphs_by_h[h]:
            capacities,edges,_area=system(graph,beta,alpha,spec)
            run=closed_baths(capacities,edges,times)
            right.append(run["u_R"])
            mass_errors.append(run["mass_relative_error"])
            minima.append(run["minimum_u"])
        mean=np.mean(np.vstack(right),axis=0)
        predictions[h]=q*mean
        diagnostics[h]={"max_mass_relative_error":max(mass_errors),"minimum_u":min(minima),
                        "u_R_seed_min":np.min(right,axis=0).tolist(),
                        "u_R_seed_max":np.max(right,axis=0).tolist()}
    return predictions, diagnostics


def residuals(curves_by_h, predictions):
    residual=[]
    for h in H_VALUES:
        points=curves_by_h[h]
        baseline=points[0][1]
        obs=np.array([k-baseline for _,k,_ in points[1:]])
        pred=np.asarray(predictions[h][1:])
        residual.extend((np.arcsinh(pred)-np.arcsinh(obs))/math.sqrt(len(obs)*len(H_VALUES)))
    return np.array(residual)


def init_forward_worker():
    global _WORKER_SPEC,_WORKER_CURVES,_WORKER_GRAPHS
    _WORKER_SPEC=json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    _WORKER_CURVES=curves()
    _WORKER_GRAPHS={(h,g["seed"]):g for h,group in graphs().items() for g in group}


def _forward_one(task):
    h,seed,alpha,beta=task
    times=np.array([point[0] for point in _WORKER_CURVES[h]])*60.0
    capacities,edges,_area=system(_WORKER_GRAPHS[h,seed],beta,alpha,_WORKER_SPEC)
    run=closed_baths(capacities,edges,times)
    return h,seed,run["u_R"],run["mass_relative_error"],run["minimum_u"]


def predict_parallel(pool: ThreadPoolExecutor, alpha: float, beta: float, q: float):
    grouped={h:[] for h in H_VALUES}
    tasks=[(h,seed,alpha,beta) for h in H_VALUES for seed in range(7001,7009)]
    for h,seed,right,mass_error,minimum_u in pool.map(_forward_one,tasks):
        grouped[h].append((seed,right,mass_error,minimum_u))
    predictions,diagnostics={},{}
    for h in H_VALUES:
        group=sorted(grouped[h],key=lambda r:r[0])
        assert [r[0] for r in group]==list(range(7001,7009))
        right=np.vstack([r[1] for r in group])
        predictions[h]=q*np.mean(right,axis=0)
        diagnostics[h]={"max_mass_relative_error":max(r[2] for r in group),
                        "minimum_u":min(r[3] for r in group),
                        "u_R_seed_min":np.min(right,axis=0).tolist(),
                        "u_R_seed_max":np.max(right,axis=0).tolist()}
    return predictions,diagnostics


def main():
    spec=json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    data=curves();nets=graphs()
    start_rows=[]
    probes=[]
    for model in ("N0","NS"):
        initial=starts(model)
        assert initial.shape==(12,2 if model=="N0" else 3)
        for index,vector in enumerate(initial):
            start_rows.append({"model":model,"start_index":index,
                               "log10_alpha":vector[0],"log10_Q":vector[1],
                               "beta":vector[2] if model=="NS" else 0.0})
        vector=initial[0]
        alpha,q=10**vector[0],10**vector[1]
        beta=vector[2] if model=="NS" else 0.0
        tick=time.perf_counter()
        predicted,diagnostics=predict(data,nets,spec,alpha,beta,q)
        r=residuals(data,predicted)
        assert len(r)==36 and np.all(np.isfinite(r))
        assert max(d["max_mass_relative_error"] for d in diagnostics.values())<=1e-8
        assert min(d["minimum_u"] for d in diagnostics.values())>=-1e-12
        probes.append({"model":model,"alpha_m_inv":alpha,"Q_original_unit":q,"beta":beta,
                       "loss_L":float(r@r),"E_sqrt_L":float(np.linalg.norm(r)),
                       "seconds":time.perf_counter()-tick,"predicted_increment_by_h":
                       {str(h):predicted[h].tolist() for h in H_VALUES},
                       "diagnostics_by_h":{str(h):diagnostics[h] for h in H_VALUES}})
        print(f"D1 {model} midpoint: L={probes[-1]['loss_L']:.6g}; "
              f"forward seconds={probes[-1]['seconds']:.3g}",flush=True)
    with (OUT/"v7_d1_starts.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(start_rows[0]));writer.writeheader();writer.writerows(start_rows)
    (OUT/"v7_d1_objective_probe.json").write_text(json.dumps({"probes":probes,
        "loss_definition":"mean over six membranes of mean over nonzero points of [asinh(pred)-asinh(obs)]^2; s0=1 μS/cm",
        "observation_role":"B Book1 only; 36 nonzero points; S geometry; eight fixed seed means"},
        ensure_ascii=False,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
