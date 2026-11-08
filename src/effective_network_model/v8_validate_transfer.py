"""Frozen v8.0.1 L conditional transfer, common mixing, and compensation rows."""

import math
from pathlib import Path

import numpy as np

from v8_compare import conductivity_array, observations, predicted_activity, read_csv, setup, write_csv

ROOT=Path(__file__).resolve().parents[2]
PACKAGE=ROOT/"原始资料/赵珈绮/16_物理诊断与接续任务_v8.0.1"
if not PACKAGE.exists():
    PACKAGE=ROOT/"taskbook"
OUT=ROOT/"WP5"


def l_observations():
    raw=read_csv(PACKAGE/"inputs/salt/salt_transport_tidy.csv")
    roles={r["curve_id"]:r["role"] for r in read_csv(PACKAGE/"inputs/salt/curve_roles_v7.csv")}
    groups={200:[],350:[],700:[]}
    for r in raw:
        if roles[r["curve_id"]]=="L":
            groups[int(r["sample_id"].split()[0])].append(r)
    assert sum(map(len,groups.values()))==70
    data={}
    for h,group in groups.items():
        ordered=sorted(group,key=lambda r:float(r["time"]))
        assert float(ordered[0]["time"])==0
        k0=float(ordered[0]["conductivity"])
        data[h]=[{"curve_id":r["curve_id"],"h_nm":h,"t_min":float(r["time"]),
                  "kappa0_uScm":k0,"delta_obs_uScm":float(r["conductivity"])-k0}
                 for r in ordered]
    return data


def main():
    new,_=setup()
    selected={(r["model"],r["split"]):r for r in read_csv(ROOT/"WP4/v8_最优参数.csv")}
    main_p={(r["model"],int(r["material_index"]),int(r["h_nm"]),int(r["seed"])):float(r["P_m_s"])
            for r in read_csv(ROOT/"WP4/v8_稳态表.csv")}
    seeds=new["geometry"]["main_seeds"]
    l_data=l_observations()
    alpha_values=new["cell"]["alpha_m_inv_grid"]
    tau_values=new["cell"]["tau_m_s_grid"]
    fit_rows=[]; chosen={}; chosen_predictions={}
    for model in ("Ncal","Scal","Tcal"):
        choice=selected[(model,"full_B")]
        local="Scal" if model=="Ncal" else model
        index=0 if model=="Ncal" else int(choice["material_index"])
        p={h:np.array([main_p[(local,index,h,seed)] for seed in seeds]) for h in l_data}
        for alpha in alpha_values:
            for tau in tau_values:
                curves={}
                errors=[]
                for h,points in l_data.items():
                    times=np.array([r["t_min"]*60 for r in points])
                    obs=np.array([r["delta_obs_uScm"] for r in points])
                    activity=predicted_activity(p[h][:,None],alpha,tau,times[None,:])
                    pred=conductivity_array(600*activity).mean(axis=0)
                    curves[h]=pred
                    mask=np.array([0<r["t_min"]<=240 for r in points])
                    if mask.any():
                        errors.append(np.mean((np.arcsinh(pred[mask])-np.arcsinh(obs[mask]))**2))
                E=math.sqrt(float(np.mean(errors)))
                fit_rows.append({"model":model,"alpha_L_m_inv":alpha,"tau_L_s":tau,
                                 "E_L_train":E,"status":"enumerated_material_locked"})
                previous=chosen.get(model)
                if previous is None or E<previous["E_L_train"]-1e-10 or (
                    abs(E-previous["E_L_train"])<=1e-10 and (tau,alpha)<(previous["tau_L_s"],previous["alpha_L_m_inv"])):
                    chosen[model]={"model":model,"material_q":choice["q"],"material_beta":choice["beta"],
                                   "alpha_L_m_inv":alpha,"tau_L_s":tau,"E_L_train":E,
                                   "time_evidence_level":"inferred_minute_unconfirmed"}
                    chosen_predictions[model]=curves
    assert len(fit_rows)==3*61*3
    write_csv(OUT/"v8_L装置网格.csv",fit_rows)
    write_csv(OUT/"v8_L选定装置参数.csv",list(chosen.values()))
    predictions=[]; summary=[]
    for model in ("Ncal","Scal","Tcal"):
        curves=chosen_predictions[model]
        late_errors=[]; terminal={}
        for h,points in l_data.items():
            pred=curves[h]
            terminal[h]=float(pred[-1]/points[-1]["delta_obs_uScm"])
            for i,point in enumerate(points):
                classification="initial" if point["t_min"]==0 else "train" if point["t_min"]<=240 else "test"
                predictions.append({**point,"model":model,"train_or_test":classification,
                                    "delta_pred_uScm":float(pred[i]),
                                    "kappa_pred_uScm":float(pred[i])+point["kappa0_uScm"],
                                    "alpha_L_m_inv":chosen[model]["alpha_L_m_inv"],
                                    "tau_L_s":chosen[model]["tau_L_s"],
                                    "time_evidence_level":"inferred_minute_unconfirmed"})
            mask=np.array([r["t_min"]>240 for r in points])
            if mask.any():
                obs=np.array([r["delta_obs_uScm"] for r in points])
                late_errors.append(np.mean((np.arcsinh(pred[mask])-np.arcsinh(obs[mask]))**2))
        summary.append({"model":model,"E_L_train":chosen[model]["E_L_train"],
                        "E_L_late":math.sqrt(float(np.mean(late_errors))),
                        "L200_terminal_pred_over_obs":terminal[200],
                        "L350_terminal_pred_over_obs":terminal[350],
                        "L700_terminal_pred_over_obs":terminal[700],
                        "time_evidence_level":"inferred_minute_unconfirmed",
                        "status":"conditional_transfer_no_material_refit"})
    assert len(predictions)==3*70
    write_csv(OUT/"v8_L逐点条件转移.csv",predictions)
    write_csv(OUT/"v8_L条件转移汇总.csv",summary)

    # Change only common mixing time at the fixed Tcal full-B material and alpha.
    b_points=observations(new["geometry"]["B_h_nm"])
    full=selected[("Tcal","full_B")]
    q_index=int(full["material_index"])
    mix=[]
    for tau in (0,300,1200):
        early=[]
        for hi,h in enumerate(new["geometry"]["B_h_nm"]):
            pts=b_points[hi*7:(hi+1)*7]
            p=np.array([main_p[("Tcal",q_index,h,seed)] for seed in seeds])
            times=np.array([r["t_min"]*60 for r in pts])
            obs=np.array([r["delta_obs_uScm"] for r in pts])
            pred=conductivity_array(600*predicted_activity(p[:,None],float(full["alpha_m_inv"]),tau,times[None,:])).mean(axis=0)
            mask=(times>0)&(times<=3600)
            residual=np.arcsinh(pred[mask])-np.arcsinh(obs[mask])
            early.append(np.mean(residual*residual))
            mix.append({"h_nm":h,"tau_s":tau,"E_early_film":math.sqrt(float(np.mean(residual*residual))),
                        "max_abs_early_residual_uScm":float(np.max(np.abs(pred[mask]-obs[mask]))),
                        "q_fixed":full["q"],"alpha_fixed_m_inv":full["alpha_m_inv"]})
    write_csv(OUT/"v8_共同混合早段.csv",mix)

    compensation=[]
    r0=0.35; q0=float(full["q"]); alpha0=float(full["alpha_m_inv"])
    for rs in (0.25,0.45):
        b=rs/r0
        mapped_q=b*q0
        mapped_alpha=alpha0/(b*b)
        compensation.append({"r_s_nm":rs,"q_mapped":mapped_q,"alpha_mapped_m_inv":mapped_alpha,
                             "within_q_box":0.003<=mapped_q<=0.5,
                             "within_alpha_box":1<=mapped_alpha<=100,
                             "scope":"approximate throat-dominated uncapped limit; port access breaks exact equivalence"})
    write_csv(OUT/"v8_组合参数近似映射.csv",compensation)
    print("L conditional fits=549; predictions=210; mixing=18; compensation=2")


if __name__=="__main__":
    main()
