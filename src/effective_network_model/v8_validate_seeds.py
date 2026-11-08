"""Frozen v8.0.1 section 5.1: independent seeds, parameters locked."""

import csv
import json
import math
from pathlib import Path

import numpy as np

from v7_network_geometry import generate, morphology
from v8_compare import conductivity_array, observations, predicted_activity, read_csv, setup, write_csv
from v8_transport import steady, weighted_graph

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP5"


def main():
    new, old = setup()
    hs = new["geometry"]["B_h_nm"]
    seeds = new["geometry"]["validation_seeds"]
    assert seeds == list(range(7009,7033))
    selected = {(r["model"],r["split"]):r for r in read_csv(ROOT/"WP4/v8_最优参数.csv")}
    points = observations(hs)
    by_h = {h:points[i*7:(i+1)*7] for i,h in enumerate(hs)}
    max_nz = max(math.ceil(h / sum(morphology(h,"S",old["geometry"]["anchors"])[:2])) for h in hs)
    steady_rows=[]
    prediction_rows=[]
    for h in hs:
        tasks=[(model,split) for model in ("Ncal","Scal","Tcal")
               for split in ("full_B","time_B",f"film_{h}")]
        for seed in seeds:
            graph, meta = generate(h,"S",seed,old,max_nz)
            info={"s_nm":meta["s_nm"],"nz":meta["nz"],"epsilon":meta["epsilon"],"nx":12,"ny":12}
            for model,split in tasks:
                choice=selected[(model,split)]
                q,beta=float(choice["q"]),float(choice["beta"])
                local="Tcal" if model=="Tcal" else "Scal"
                edges,area=weighted_graph(graph,info,q,local,beta)
                solved=steady(len(graph["nodes"]),edges,area)
                steady_rows.append({"model":model,"split":split,"h_nm":h,"seed":seed,
                                    "q":q,"beta":beta,"G_m3s":solved["G_m3s"],
                                    "P_m_s":solved["P_m_s"],
                                    "conservation_relative":solved["conservation_relative"],
                                    "status":"PASS_TECHNICAL"})
    assert len(steady_rows)==1296
    write_csv(OUT/"v8_独立种子稳态.csv",steady_rows)
    grouped={}
    for row in steady_rows:
        grouped.setdefault((row["model"],row["split"],row["h_nm"]),[]).append(float(row["P_m_s"]))
    for (model,split,h),p_list in grouped.items():
        assert len(p_list)==24
        choice=selected[(model,split)]
        times=np.array([p["t_min"]*60 for p in by_h[h]])
        activity=predicted_activity(np.array(p_list)[:,None],float(choice["alpha_m_inv"]),
                                    float(choice["tau_s"]),times[None,:])
        each=conductivity_array(600*activity)
        mean=each.mean(axis=0)
        for index,point in enumerate(by_h[h]):
            prediction_rows.append({"model":model,"split":split,"curve_id":point["curve_id"],
                                    "h_nm":h,"t_min":point["t_min"],
                                    "delta_obs_uScm":point["delta_obs_uScm"],
                                    "delta_independent_mean_uScm":float(mean[index]),
                                    "delta_independent_min_uScm":float(each[:,index].min()),
                                    "delta_independent_max_uScm":float(each[:,index].max()),
                                    "seed_set":"7009-7032","q":choice["q"],"beta":choice["beta"],
                                    "alpha_m_inv":choice["alpha_m_inv"],"tau_s":choice["tau_s"],
                                    "status":"fixed_parameters_no_reselection"})
    assert len(prediction_rows)==378
    write_csv(OUT/"v8_独立种子预测.csv",prediction_rows)
    by={(r["model"],r["split"],int(r["h_nm"]),float(r["t_min"])):r for r in prediction_rows}
    main={(r["model"],r["split"],int(r["h_nm"]),float(r["t_min"])):r
          for r in read_csv(ROOT/"WP4/v8_全部逐点预测.csv")}
    summary=[]
    for model in ("Ncal","Scal","Tcal"):
        full_errors=[]; time_errors=[]; film_errors=[]; ratio=[]; shifts=[]
        for h in hs:
            for point in by_h[h][1:]:
                t=point["t_min"]
                obs=point["delta_obs_uScm"]
                for split,collector in (("full_B",full_errors),("time_B",time_errors),
                                        (f"film_{h}",film_errors)):
                    if split=="time_B" and t<=60:
                        continue
                    rec=by[(model,split,h,t)]
                    pred=float(rec["delta_independent_mean_uScm"])
                    collector.append((h,(math.asinh(pred)-math.asinh(obs))**2))
                    primary=float(main[(model,split,h,t)]["delta_pred_uScm"])
                    shifts.append(abs(pred-primary))
            ratio.append(float(by[(model,f"film_{h}",h,by_h[h][-1]["t_min"])]["delta_independent_mean_uScm"])/by_h[h][-1]["delta_obs_uScm"])
        def equal_film(items):
            return math.sqrt(np.mean([np.mean([value for hh,value in items if hh==h]) for h in hs]))
        summary.append({"model":model,"E_fullB_independent":equal_film(full_errors),
                        "E_time_independent":equal_film(time_errors),
                        "E_film_independent":equal_film(film_errors),
                        "film_2600_endpoint_pred_over_obs":ratio[-1],
                        "film_endpoint_ratio_min":min(ratio),"film_endpoint_ratio_max":max(ratio),
                        "max_abs_primary_to_independent_uScm":max(shifts),
                        "seed_count":24,"linear_solves":432,
                        "interpretation":"geometry seed variability, not experimental uncertainty"})
    write_csv(OUT/"v8_独立种子汇总.csv",summary)
    print(f"independent seed solves={len(steady_rows)}; predictions={len(prediction_rows)}")


if __name__=="__main__":
    main()
