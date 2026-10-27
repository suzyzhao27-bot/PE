"""Apply frozen v8.0.1 science gates to saved results; no new fit."""

import json
import math
from pathlib import Path

import numpy as np

from v8_compare import conductivity_array, observations, predicted_activity, read_csv, setup

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"WP5"


def main():
    new,_=setup()
    models={r["model"]:r for r in read_csv(ROOT/"WP4/v8_模型比较.csv")}
    technical=read_csv(ROOT/"WP3/v8_技术检查.csv")
    selected={(r["model"],r["split"]):r for r in read_csv(ROOT/"WP4/v8_最优参数.csv")}
    t=models["Tcal"]; n=models["Ncal"]; s=models["Scal"]
    et,ef=float(t["E_time"]),float(t["E_film"])
    improvements={"time_vs_Ncal":1-et/float(n["E_time"]),
                  "film_vs_Ncal":1-ef/float(n["E_film"]),
                  "time_vs_Scal":1-et/float(s["E_time"]),
                  "film_vs_Scal":1-ef/float(s["E_film"])}
    seed={r["model"]:r for r in read_csv(OUT/"v8_独立种子汇总.csv")}
    domain=read_csv(OUT/"v8_域尺寸汇总.csv")
    t_domain=[r for r in domain if r["model"]=="Tcal"]
    scenarios=read_csv(OUT/"v8_固定十情景汇总.csv")
    water=read_csv(OUT/"v8_水力上界逐种子.csv")
    transfer={r["model"]:r for r in read_csv(OUT/"v8_L条件转移汇总.csv")}
    rows=read_csv(ROOT/"WP4/v8_全网格损失.csv")
    full=selected[("Tcal","full_B")]
    min_e=float(full["E_train"])
    near=[r for r in rows if r["model"]=="Tcal" and r["split"]=="full_B" and
          r["status"]!="NOT_STARTED_BY_GATE:G2" and float(r["E_train"])<=1.05*min_e+0.01]
    near_q=sorted(set(float(r["q"]) for r in near))
    q_grid=new["models"]["Tcal"]["q_grid"]
    q0=float(full["q"])
    selected_index=next(i for i,q in enumerate(q_grid) if math.isclose(q,q0,rel_tol=1e-14))
    adjacent=[q_grid[i] for i in (selected_index-1,selected_index+1) if 0<=i<len(q_grid)]
    hs=new["geometry"]["B_h_nm"]
    seeds=new["geometry"]["main_seeds"]
    steady={(float(r["q"]),int(r["h_nm"]),int(r["seed"])):float(r["P_m_s"])
            for r in read_csv(ROOT/"WP4/v8_稳态表.csv") if r["model"]=="Tcal"}
    points=observations(hs)
    obs=np.array([[r["delta_obs_uScm"] for r in points[i*7:(i+1)*7]] for i in range(6)])
    times=np.array([[r["t_min"]*60 for r in points[i*7:(i+1)*7]] for i in range(6)])
    adjacent_report=[]
    for q in [q0]+adjacent:
        if not all((q,h,seed) in steady for h in hs for seed in seeds):
            adjacent_report.append({"q":q,"status":"NOT_STARTED_BY_GATE:G2"})
            continue
        p=np.array([[steady[(q,h,seed)] for seed in seeds] for h in hs])
        activity=predicted_activity(p[:,:,None],float(full["alpha_m_inv"]),float(full["tau_s"]),times[:,None,:])
        pred=conductivity_array(600*activity).mean(axis=1)
        residual=np.arcsinh(pred[:,1:])-np.arcsinh(obs[:,1:])
        error=math.sqrt(float(np.mean(np.mean(residual*residual,axis=1))))
        adjacent_report.append({"q":q,"E_fullB_fixed_cell":error,
                                "endpoint_2600_pred_over_obs":float(pred[-1,-1]/obs[-1,-1]),
                                "endpoint_200_over_350_pred":float(pred[2,-1]/pred[3,-1]),
                                "status":"computed_from_saved_steady"})
    endpoint_ratios=[float(r["endpoint_pred_over_obs"]) for r in read_csv(ROOT/"WP4/v8_逐膜残差.csv")
                     if r["model"]=="Tcal" and r["split"].startswith("film_") and
                     int(r["h_nm"])==int(r["split"].split("_")[1])]
    prediction_rows=read_csv(ROOT/"WP4/v8_全部逐点预测.csv")
    def endpoint_ratio_for(split):
        endpoint={int(r["h_nm"]):float(r["delta_pred_uScm"]) for r in prediction_rows
                  if r["model"]=="Tcal" and r["split"]==split and
                  math.isclose(float(r["t_min"]),120) and int(r["h_nm"]) in (200,350)}
        return endpoint[200]/endpoint[350]
    result={
        "scope":"v8.0.1 frozen retrospective validation; no new fit",
        "technical_8_of_8":len(technical)==8 and all(r["status"]=="PASS_TECHNICAL" for r in technical),
        "relative_improvements":improvements,
        "science_gates":{
            "Ncal_20pct_both":improvements["time_vs_Ncal"]>=0.2 and improvements["film_vs_Ncal"]>=0.2,
            "Scal_10pct_both":improvements["time_vs_Scal"]>=0.1 and improvements["film_vs_Scal"]>=0.1,
            "E_time_and_film_at_most_0p70":et<=0.7 and ef<=0.7,
            "all_film_endpoints_1over3_to3":all(1/3<=x<=3 for x in endpoint_ratios),
            "Tcal_domain_all_within_10pct":all(r["within_10pct"]=="True" for r in t_domain),
            "conditional_water_bound_all_seeds_below_5pct_120min":all(float(r["volume_fraction_upper_120min"])<0.05 for r in water),
        },
        "Tcal_E_time":et,"Tcal_E_film":ef,
        "Tcal_film_2600_endpoint_pred_over_obs":float(t["film_2600_endpoint_pred_over_obs"]),
        "Tcal_200_over_350_fullB":float(t["fullB_200_over_350_pred"]),
        "Tcal_200_over_350_film200":endpoint_ratio_for("film_200"),
        "Tcal_200_over_350_film350":endpoint_ratio_for("film_350"),
        "observed_200_over_350":float(t["fullB_200_over_350_obs"]),
        "Tcal_independent_E_time":float(seed["Tcal"]["E_time_independent"]),
        "Tcal_independent_E_film":float(seed["Tcal"]["E_film_independent"]),
        "Tcal_domain_pass_count":sum(r["within_10pct"]=="True" for r in t_domain),
        "Tcal_domain_total":len(t_domain),
        "Scal_domain_exception":"full_B 350 nm, 13.5685% mean P change",
        "water_max_seed_fraction_120min":max(float(r["volume_fraction_upper_120min"]) for r in water),
        "fixed_scenario_2600_ratio_range":[min(float(r["film_2600_endpoint_pred_over_obs"]) for r in scenarios),
                                            max(float(r["film_2600_endpoint_pred_over_obs"]) for r in scenarios)],
        "fixed_scenario_chi_2600_350_range":[min(float(r["chi_2600_over_350"]) for r in scenarios),
                                               max(float(r["chi_2600_over_350"]) for r in scenarios)],
        "conditional_L_Tcal_E_late":float(transfer["Tcal"]["E_L_late"]),
        "near_optimal_rule":"E_train <= 1.05 E_min + 0.01; not a confidence set",
        "near_optimal_grid_count":len(near),"near_optimal_q_values":near_q,
        "adjacent_q_fixed_cell":adjacent_report,
        "decision":"SCIENTIFIC_NOT_SUPPORTED",
        "interpretation":"Failure of frozen joint prediction; does not exclude all pore physics or validate a different mechanism.",
    }
    (OUT/"v8_主张判定.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(f"science gates={sum(result['science_gates'].values())}/{len(result['science_gates'])}; decision={result['decision']}")


if __name__=="__main__":
    main()
