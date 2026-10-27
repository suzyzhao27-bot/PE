"""Report one completed 12-start full-B fit without converting train to validation."""

from __future__ import annotations

import argparse
import csv
import json
import math

import numpy as np

from v7_d1_objective import H_VALUES, OUT, curves, graphs, predict, residuals
from v7_network_transport import SPEC_PATH


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("model",choices=("N0","NS"))
    args=parser.parse_args()
    records=[json.loads(line) for line in (OUT/"v7_d1_fit_fullB_runs.jsonl").read_text(encoding="utf-8").splitlines()]
    selected=[r for r in records if r["model"]==args.model and r["task"]=="full_B"]
    if len(selected)!=12 or sorted(r["start_index"] for r in selected)!=list(range(12)):
        raise RuntimeError("12 unique frozen starts have not finished")
    good=[r for r in selected if r["status"]=="solver_success" and r["best_L"] is not None]
    stable=any(abs(a["best_L"]-b["best_L"])/min(a["best_L"],b["best_L"])<=0.01
               for i,a in enumerate(good) for b in good[i+1:])
    best=min(good,key=lambda r:r["best_L"])
    x=best["best_transformed"]
    alpha,q=10**x[0],10**x[1]
    beta=x[2] if args.model=="NS" else 0.0
    spec=json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    data=curves();nets=graphs()
    predicted,diagnostics=predict(data,nets,spec,alpha,beta,q)
    residual=residuals(data,predicted)
    recomputed=float(residual@residual)
    assert abs(recomputed-best["best_L"])<=1e-6*max(1.0,best["best_L"])
    rows=[]
    per_membrane=[]
    for h in H_VALUES:
        points=data[h]
        baseline=points[0][1]
        obs=[];pred=[];rs=[]
        for index,((t,k,curve_id),y) in enumerate(zip(points,predicted[h],strict=True)):
            measured=k-baseline
            r=math.asinh(y)-math.asinh(measured)
            rows.append({"model":args.model,"series":"B","h_nm":h,"curve_id":curve_id,
                         "point_index":index,"time_min":t,"conductivity_uS_cm":k,
                         "baseline_uS_cm":baseline,"observed_increment_uS_cm":measured,
                         "predicted_increment_uS_cm":y,
                         "prediction_seed_min_uS_cm":q*diagnostics[h]["u_R_seed_min"][index],
                         "prediction_seed_max_uS_cm":q*diagnostics[h]["u_R_seed_max"][index],
                         "asinh_residual":r,"in_training":index>0})
            if index>0:
                obs.append(measured);pred.append(y);rs.append(r)
        per_membrane.append({"h_nm":h,"n_train":len(obs),
                             "RMSE_original_uS_cm":float(np.sqrt(np.mean((np.array(pred)-np.array(obs))**2))),
                             "E_asinh":float(np.sqrt(np.mean(np.array(rs)**2))),
                             "endpoint_ratio_pred_over_obs":pred[-1]/obs[-1] if obs[-1]!=0 else None,
                             "max_mass_relative_error":diagnostics[h]["max_mass_relative_error"],
                             "minimum_u":diagnostics[h]["minimum_u"]})
    with (OUT/f"v7_{args.model}_fullB_points.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    p200=next(row for row in per_membrane if row["h_nm"]==200)
    p350=next(row for row in per_membrane if row["h_nm"]==350)
    endpoint_obs={h:next(row for row in rows if row["h_nm"]==h and row["point_index"]==6)["observed_increment_uS_cm"]
                  for h in (200,350)}
    endpoint_pred={h:next(row for row in rows if row["h_nm"]==h and row["point_index"]==6)["predicted_increment_uS_cm"]
                   for h in (200,350)}
    bounds=[(0.0,5.0),(0.0,6.0)]+([(0.0,6.0)] if args.model=="NS" else [])
    boundary=[name for name,value,(lo,hi) in zip(("log10_alpha","log10_Q","beta"),x,bounds)
              if min(value-lo,hi-value)<0.01*(hi-lo)]
    summary={"model":args.model,"status":"multi_start_stable" if stable else "local_optimum_unstable",
             "12_start_status_counts":{status:sum(r["status"]==status for r in selected)
                                       for status in sorted(set(r["status"] for r in selected))},
             "actual_objective_calls_total":sum(r["actual_objective_calls"] for r in selected),
             "best_start_index":best["start_index"],"best_transformed":x,
             "alpha_m_inv":alpha,"Q_uS_cm_per_unit_uR":q,"beta":beta,
             "boundary_within_one_percent":boundary,
             "L_recomputed":recomputed,"E_sqrt_L":math.sqrt(recomputed),
             "best_vs_worst_success_L_relative":
             (max(r["best_L"] for r in good)-min(r["best_L"] for r in good))/min(r["best_L"] for r in good),
             "per_membrane":per_membrane,
             "worst_membrane_by_asinh_E":max(per_membrane,key=lambda r:r["E_asinh"])["h_nm"],
             "endpoint_350_over_200_observed":endpoint_obs[350]/endpoint_obs[200],
             "endpoint_350_over_200_predicted":endpoint_pred[350]/endpoint_pred[200],
             "scope":"full B training only; no time or whole-membrane holdout"}
    (OUT/f"v7_{args.model}_fullB_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(f"{args.model} full B: stable={stable}; E={summary['E_sqrt_L']:.6g}; "
          f"worst h={summary['worst_membrane_by_asinh_E']} nm; "
          f"350/200 predicted={summary['endpoint_350_over_200_predicted']:.6g}")


if __name__=="__main__":
    main()
