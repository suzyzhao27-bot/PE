"""Frozen v8.0.1 steady table and equal-film B-grid comparison."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from v7_network_geometry import generate, morphology
from v8_transport import baths, steady, weighted_graph

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/16_物理诊断与接续任务_v8.0.1"
if not PACKAGE.exists():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP4"
STEADY = OUT / "v8_稳态表.csv"
SPLITS = ("full_B", "time_B", "film_65", "film_120", "film_200", "film_350", "film_700", "film_2600")


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def setup():
    new = json.loads((PACKAGE / "configs/frozen_spec_v8.json").read_text(encoding="utf-8"))
    old = json.loads((PACKAGE / "inputs/legacy_v7_spec.json").read_text(encoding="utf-8"))
    old["geometry"] = dict(old["geometry"])
    old["geometry"].update(nx=12, ny=12, p_z=0.65, sigma_ln=0.35, q=0.5)
    return new, old


def table():
    new, old = setup()
    pre = json.loads((ROOT / "WP2/预筛可行集.json").read_text(encoding="utf-8"))
    hs = new["geometry"]["B_h_nm"]
    seeds = new["geometry"]["main_seeds"]
    q_values = pre["accepted_q"]
    beta_values = new["models"]["Scal"]["beta_grid"]
    max_nz = max(math.ceil(h / sum(morphology(h,"S",old["geometry"]["anchors"])[:2])) for h in hs)
    rows = []
    for h in hs:
        for seed in seeds:
            graph, meta = generate(h,"S",seed,old,max_nz)
            info = {"s_nm":meta["s_nm"],"nz":meta["nz"],"epsilon":meta["epsilon"],"nx":12,"ny":12}
            for model, values in (("Scal", beta_values), ("Tcal", q_values)):
                for material_index, value in enumerate(values):
                    q = 0.5 if model == "Scal" else value
                    beta = value if model == "Scal" else 0
                    edges, area = weighted_graph(graph,info,q,model,beta)
                    solved = steady(len(graph["nodes"]),edges,area)
                    assert solved["G_m3s"] >= 0
                    assert solved["conservation_relative"] < 1e-7
                    rows.append({"model":model,"material_index":material_index,
                                 "beta":beta,"q":q,"h_nm":h,"seed":seed,
                                 "G_m3s":solved["G_m3s"],"P_m_s":solved["P_m_s"],
                                 "ion_edges":len(edges),"spanning_nodes":solved["spanning_nodes"],
                                 "linear_residual_relative":solved["residual_relative"],
                                 "conservation_relative":solved["conservation_relative"],
                                 "status":"PASS_TECHNICAL"})
    assert len(rows) == (len(beta_values)+len(q_values))*6*8
    write_csv(STEADY,rows)
    print(f"steady solves={len(rows)}; Scal={len(beta_values)*48}; Tcal={len(q_values)*48}")


def predicted_activity(p, alpha, tau, times):
    """Vectorized exact receiving-bulk solution for frozen equal-volume cell."""
    r = alpha * p
    if tau == 0:
        return -np.expm1(-2*r*times)/2
    b = r*(1+1/0.05) + 1/tau
    c = 2*r/tau
    root = np.sqrt(np.maximum(b*b-4*c,0))
    fast = -(b+root)/2
    slow = np.divide(-2*c,b+root,out=np.zeros_like(b),where=b+root>0)
    denominator = 2*(slow-fast)
    return 0.5 + (fast*np.exp(slow*times)-slow*np.exp(fast*times))/denominator


def conductivity_array(c):
    join = 10**(0.0165+0.9398-3.2793)*1000/58.44
    c = np.maximum(c,0)
    result = np.zeros_like(c)
    low = (c>0)&(c<join)
    result[low] = 10*c[low]/join
    high = c>=join
    y = np.log10(c[high]*58.44/1000)
    discriminant = 0.9398**2 - 4*0.0165*(-3.2793-y)
    result[high] = 10**((-0.9398+np.sqrt(discriminant))/(2*0.0165))
    return result


def observations(hs):
    raw = read_csv(PACKAGE / "inputs/salt/salt_transport_tidy.csv")
    role = {r["curve_id"]:r["role"] for r in read_csv(PACKAGE / "inputs/salt/curve_roles_v7.csv")}
    b = [r for r in raw if role[r["curve_id"]] == "B"]
    by_h = {h:[] for h in hs}
    for row in b:
        h = int(row["sample_id"].split()[0])
        by_h[h].append(row)
    rows = []
    for h in hs:
        group = sorted(by_h[h], key=lambda r:float(r["time"]))
        assert len(group)==7 and float(group[0]["time"])==0
        k0 = float(group[0]["conductivity"])
        for row in group:
            rows.append({"curve_id":row["curve_id"],"h_nm":h,
                         "t_min":float(row["time"]),"kappa0_uScm":k0,
                         "kappa_obs_uScm":float(row["conductivity"]),
                         "delta_obs_uScm":float(row["conductivity"])-k0})
    assert len(rows)==42 and sum(r["t_min"]>0 for r in rows)==36
    return rows


def training_mask(split, hs, times):
    nonzero = times>0
    if split=="time_B":
        return nonzero&(times<=60)
    if split.startswith("film_"):
        held=int(split.split("_")[1])
        return nonzero&(np.array(hs)[:,None]!=held)
    return nonzero


def equal_film_error(predicted, observed, mask):
    residual = np.arcsinh(predicted)-np.arcsinh(observed)
    return float(np.sqrt(np.mean([np.mean(residual[h,mask[h]]**2)
                                  for h in range(mask.shape[0]) if np.any(mask[h])])))


def selection_key(row, model):
    value = row["beta"] if model in ("Ncal","Scal") else -row["q"]
    return (row["tau_s"],value,row["alpha_m_inv"])


def compare():
    from sys import path
    path.insert(0, str(ROOT))
    from src.data_ingest.v8_step1_units import concentration, segment
    for p_probe in (1e-12, 1e-8, 1e-4):
        for alpha_probe in (1, 100):
            for tau_probe in (0, 300, 1200):
                t_probe = np.array([0, 1200, 3600, 7200], dtype=float)
                fast = predicted_activity(np.array(p_probe), alpha_probe, tau_probe, t_probe)
                reference = baths(p_probe, alpha_probe, tau_probe, t_probe)["receiver_c_molm3"]/600
                assert np.max(np.abs(fast-reference)) < 1e-9
    new,_ = setup()
    hs = new["geometry"]["B_h_nm"]
    seeds = new["geometry"]["main_seeds"]
    source = read_csv(STEADY)
    steady = {}
    for r in source:
        key=(r["model"],int(r["material_index"]),int(r["h_nm"]),int(r["seed"]))
        steady[key]=float(r["P_m_s"])
    beta_values = new["models"]["Scal"]["beta_grid"]
    prescreen = json.loads((ROOT / "WP2/预筛可行集.json").read_text(encoding="utf-8"))
    q_values = prescreen["accepted_q"]
    excluded_q = [q for q in new["models"]["Tcal"]["q_grid"] if q not in q_values]
    points = observations(hs)
    times = np.array([[r["t_min"]*60 for r in points[i*7:(i+1)*7]] for i in range(6)])
    observed = np.array([[r["delta_obs_uScm"] for r in points[i*7:(i+1)*7]] for i in range(6)])
    masks={split:training_mask(split,hs,times/60) for split in SPLITS}
    alpha_values=new["cell"]["alpha_m_inv_grid"]
    tau_values=new["cell"]["tau_m_s_grid"]
    fit_rows=[]
    selected={}
    predictions={}
    for model, source_model, values in (("Ncal","Scal",[0]),("Scal","Scal",beta_values),
                                        ("Tcal","Tcal",q_values)):
        for material_index,value in enumerate(values):
            index = 0 if model=="Ncal" else material_index
            p = np.array([[steady[(source_model,index,h,seed)] for seed in seeds] for h in hs])
            for alpha in alpha_values:
                for tau in tau_values:
                    u = predicted_activity(p[:,:,None],alpha,tau,times[:,None,:])
                    delta = np.mean(conductivity_array(600*u),axis=1)
                    assert np.all(np.isfinite(delta)) and np.min(delta)>=-1e-8
                    for split in SPLITS:
                        error=equal_film_error(delta,observed,masks[split])
                        row={"model":model,"split":split,"material_index":index,
                             "q":0.5 if model!="Tcal" else value,
                             "beta":value if model=="Scal" else 0,
                             "alpha_m_inv":alpha,"tau_s":tau,"E_train":error,
                             "status":"enumerated"}
                        fit_rows.append(row)
                        old=selected.get((model,split))
                        if old is None or error<old["E_train"]-1e-10 or (
                            abs(error-old["E_train"])<=1e-10 and selection_key(row,model)<selection_key(old,model)):
                            selected[(model,split)]=row
                            predictions[(model,split)]=delta.copy()
    assert len(selected)==24
    selected_ids={(r["model"],r["split"],r["material_index"],r["alpha_m_inv"],r["tau_s"])
                  for r in selected.values()}
    for row in fit_rows:
        if (row["model"],row["split"],row["material_index"],row["alpha_m_inv"],row["tau_s"]) in selected_ids:
            row["status"]="selected"
    for q in excluded_q:
        for alpha in alpha_values:
            for tau in tau_values:
                for split in SPLITS:
                    fit_rows.append({"model":"Tcal","split":split,"material_index":"",
                                     "q":q,"beta":0,"alpha_m_inv":alpha,"tau_s":tau,
                                     "E_train":"","status":"NOT_STARTED_BY_GATE:G2"})
    write_csv(OUT/"v8_全网格损失.csv",fit_rows)
    write_csv(OUT/"v8_最优参数.csv",list(selected.values()))
    prediction_rows=[]
    residual_rows=[]
    spec_sha256 = hashlib.sha256((PACKAGE/"configs/frozen_spec_v8.json").read_bytes()).hexdigest()
    for (model,split),row in selected.items():
        delta=predictions[(model,split)]
        mask=masks[split]
        for hi,h in enumerate(hs):
            per=[]
            for ti,point in enumerate(points[hi*7:(hi+1)*7]):
                x=float(delta[hi,ti])
                residual=float(np.arcsinh(x)-np.arcsinh(point["delta_obs_uScm"]))
                train = bool(mask[hi,ti])
                classification="initial" if ti==0 else "train" if train else "test"
                prediction_rows.append({**point,"model":model,"split":split,
                                        "train_or_test":classification,
                                        "delta_pred_uScm":x,"kappa_pred_uScm":x+point["kappa0_uScm"],
                                        "c_conditional_molm3":concentration(point["delta_obs_uScm"]),
                                        "conversion_segment":segment(point["delta_obs_uScm"]),
                                        "prediction_conversion_segment":segment(x),
                                        "time_evidence_level":"explicit_origin_units",
                                        "asinh_residual":residual,
                                        "q":row["q"],"beta":row["beta"],
                                        "alpha_m_inv":row["alpha_m_inv"],"tau_s":row["tau_s"],
                                        "C0_molm3":600,"gamma":1,"seed_set":"7001-7008",
                                        "spec_sha256":spec_sha256})
                if ti>0:
                    per.append((train,residual,x-point["delta_obs_uScm"]))
            residual_rows.append({"model":model,"split":split,"h_nm":h,
                                  "E_train_film":math.sqrt(np.mean([r*r for train,r,_ in per if train])) if any(train for train,_,_ in per) else "",
                                  "E_test_film":math.sqrt(np.mean([r*r for train,r,_ in per if not train])) if any(not train for train,_,_ in per) else "",
                                  "MAE_uScm":np.mean([abs(diff) for _,_,diff in per]),
                                  "endpoint_pred_over_obs":float(delta[hi,-1]/observed[hi,-1])})
    write_csv(OUT/"v8_全部逐点预测.csv",prediction_rows)
    write_csv(OUT/"v8_逐膜残差.csv",residual_rows)
    comparisons=[]
    for model in ("Ncal","Scal","Tcal"):
        full=predictions[(model,"full_B")]
        time_pred=predictions[(model,"time_B")]
        time_mask=times>60*60
        fold_errors=[]
        fold_end_ratios=[]
        for hi,h in enumerate(hs):
            fold=predictions[(model,f"film_{h}")]
            residual=np.arcsinh(fold[hi,1:])-np.arcsinh(observed[hi,1:])
            fold_errors.append(float(np.mean(residual**2)))
            fold_end_ratios.append(float(fold[hi,-1]/observed[hi,-1]))
        comparisons.append({"model":model,"E_fullB":selected[(model,"full_B")]["E_train"],
                            "E_time":equal_film_error(time_pred,observed,time_mask),
                            "E_film":math.sqrt(np.mean(fold_errors)),
                            "fullB_MAE_uScm":float(np.mean(np.abs(full[:,1:]-observed[:,1:]))),
                            "fullB_200_over_350_pred":float(full[2,-1]/full[3,-1]),
                            "fullB_200_over_350_obs":float(observed[2,-1]/observed[3,-1]),
                            "fullB_2600_endpoint_pred_over_obs":float(full[-1,-1]/observed[-1,-1]),
                            "film_2600_endpoint_pred_over_obs":fold_end_ratios[-1],
                            "film_endpoint_ratio_min":min(fold_end_ratios),
                            "film_endpoint_ratio_max":max(fold_end_ratios)})
    write_csv(OUT/"v8_模型比较.csv",comparisons)
    print(f"fit grid={len(fit_rows)}; selected={len(selected)}; prediction rows={len(prediction_rows)}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=("tables","compare"))
    args=parser.parse_args()
    table() if args.stage=="tables" else compare()
