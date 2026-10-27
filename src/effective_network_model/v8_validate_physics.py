"""Frozen v8.0.1 fixed scenarios and conditional hydraulic volume bounds."""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from v7_network_geometry import generate, morphology
from v8_compare import conductivity_array, observations, predicted_activity, read_csv, setup, write_csv
from v8_transport import steady, weighted_graph

ROOT=Path(__file__).resolve().parents[2]
PACKAGE=ROOT/"原始资料/赵珈绮/16_物理诊断与接续任务_v8.0.1"
if not PACKAGE.exists():
    PACKAGE=ROOT/"taskbook"
OUT=ROOT/"WP5"


def max_nz_for(hs,geom):
    return max(math.ceil(h/sum(morphology(h,"S",geom["anchors"])[:2])) for h in hs)


def geometry_info(meta,nxy):
    return {"s_nm":meta["s_nm"],"nz":meta["nz"],"epsilon":meta["epsilon"],"nx":nxy,"ny":nxy}


def fixed_scenarios():
    new,base_spec=setup()
    hs=new["geometry"]["B_h_nm"]
    seeds=new["geometry"]["main_seeds"]
    choice=next(r for r in read_csv(ROOT/"WP4/v8_最优参数.csv")
                if r["model"]=="Tcal" and r["split"]=="full_B")
    q,alpha,tau=map(float,(choice["q"],choice["alpha_m_inv"],choice["tau_s"]))
    scenarios=[("r_s_0p25",{"r_s":0.25e-9}),("r_s_0p45",{"r_s":0.45e-9}),
               ("dp_0p8",{"dp_factor":0.8}),("dp_1p2",{"dp_factor":1.2}),
               ("sigma_0p25",{"sigma_ln":0.25}),("sigma_0p45",{"sigma_ln":0.45}),
               ("pz_0p5",{"p_z":0.5}),("pz_0p8",{"p_z":0.8}),
               ("gamma_0p8",{"gamma":0.8}),("gamma_1p2",{"gamma":1.2})]
    base_rows=read_csv(ROOT/"WP4/v8_稳态表.csv")
    base_p={(int(r["h_nm"]),int(r["seed"])):float(r["P_m_s"]) for r in base_rows
            if r["model"]=="Tcal" and int(r["material_index"])==int(choice["material_index"])}
    rows=[]
    points=observations(hs)
    obs=np.array([[r["delta_obs_uScm"] for r in points[i*7:(i+1)*7]] for i in range(6)])
    times=np.array([[r["t_min"]*60 for r in points[i*7:(i+1)*7]] for i in range(6)])
    summaries=[]
    for label,change in scenarios:
        geom=base_spec["geometry"].copy()
        if "dp_factor" in change:
            geom["anchors"]=[{**a,"dp_nm":a["dp_nm"]*change["dp_factor"]} for a in geom["anchors"]]
        if "sigma_ln" in change: geom["sigma_ln"]=change["sigma_ln"]
        if "p_z" in change: geom["p_z"]=change["p_z"]
        spec={**base_spec,"geometry":geom}
        max_nz=max_nz_for(hs,geom)
        grid=np.empty((6,8))
        for hi,h in enumerate(hs):
            for si,seed in enumerate(seeds):
                if "gamma" in change:
                    p=base_p[(h,seed)]
                    conservation="reused main steady table"
                    solve_count=0
                else:
                    graph,meta=generate(h,"S",seed,spec,max_nz)
                    info=geometry_info(meta,12)
                    edges,area=weighted_graph(graph,info,q,"Tcal",r_s=change.get("r_s",0.35e-9))
                    solved=steady(len(graph["nodes"]),edges,area)
                    p=solved["P_m_s"]
                    conservation=solved["conservation_relative"]
                    solve_count=1
                grid[hi,si]=p
                rows.append({"scenario":label,"h_nm":h,"seed":seed,"q_fixed":q,
                             "alpha_fixed_m_inv":alpha,"tau_fixed_s":tau,
                             "P_m_s":p,"conservation_relative":conservation,
                             "new_steady_solve":solve_count})
        activity=predicted_activity(grid[:,:,None],alpha,tau,times[:,None,:])
        gamma=change.get("gamma",1)
        delta=np.mean(gamma*conductivity_array(600*activity),axis=1)
        residual=np.arcsinh(delta[:,1:])-np.arcsinh(obs[:,1:])
        error=math.sqrt(float(np.mean(np.mean(residual*residual,axis=1))))
        mean_p=np.mean(grid,axis=1)
        summaries.append({"scenario":label,"E_fullB_fixed":error,
                          "film_2600_endpoint_pred_over_obs":float(delta[-1,-1]/obs[-1,-1]),
                          "film_200_over_350_end":float(delta[2,-1]/delta[3,-1]),
                          "P_2600_over_350":float(mean_p[-1]/mean_p[3]) if mean_p[3]>0 else "",
                          "chi_2600_over_350":float(mean_p[-1]*2600/(mean_p[3]*350)) if mean_p[3]>0 else "",
                          "linear_solves":sum(r["new_steady_solve"] for r in rows if r["scenario"]==label),
                          "status":"fixed_no_refit"})
    assert len(rows)==480 and sum(r["new_steady_solve"] for r in rows)==384
    write_csv(OUT/"v8_固定十情景稳态.csv",rows)
    write_csv(OUT/"v8_固定十情景汇总.csv",summaries)
    print(f"fixed scenarios={len(summaries)}; steady solves=384; rows={len(rows)}")


def hydraulic_edges(graph,meta,q):
    nodes=graph["nodes"]
    n=len(nodes)
    s=meta["s_nm"]*1e-9
    dz=graph["h_nm"]*1e-9/meta["nz"]
    epsilon=meta["epsilon"]
    nxy=12*12
    result=[]
    for e in graph["edges"]:
        i,j=e["i"],e["j"]
        if e["port_flag"]:
            node=j if i=="L" else i
            raw=q*nodes[node]["r_m"]
            face=s*s
        else:
            raw=q*min(nodes[i]["r_m"],nodes[j]["r_m"])
            face=s*s if i//nxy!=j//nxy else s*dz
        a=min(raw,math.sqrt(epsilon*face/math.pi))
        g=math.pi*a**4/(8*1e-3*e["L_m"])
        ii=n if i=="L" else n+1 if i=="R" else i
        jj=n if j=="L" else n+1 if j=="R" else j
        result.append((ii,jj,g))
    return result,(12*s)**2


def water_bounds():
    new,old=setup()
    hs=new["geometry"]["B_h_nm"]
    seeds=new["geometry"]["main_seeds"]
    choices={(r["model"],r["split"]):r for r in read_csv(ROOT/"WP4/v8_最优参数.csv")}
    max_nz=max_nz_for(hs,old["geometry"])
    rows=[]
    for h in hs:
        tasks=["full_B","time_B",f"film_{h}"]
        for seed in seeds:
            graph,meta=generate(h,"S",seed,old,max_nz)
            for split in tasks:
                choice=choices[("Tcal",split)]
                q=float(choice["q"])
                edges,area=hydraulic_edges(graph,meta,q)
                result=steady(len(graph["nodes"]),edges,area)
                H=result["G_m3s"]
                alpha=float(choice["alpha_m_inv"])
                osmotic=2*8.314462618*298.15*600
                bound120=(H/area)*osmotic*alpha*7200
                rows.append({"split":split,"h_nm":h,"seed":seed,"q":q,"alpha_m_inv":alpha,
                             "H_m3_Pa_s":H,"area_box_m2":area,"delta_pi_Pa":osmotic,
                             "volume_fraction_upper_120min":bound120,
                             "conservation_relative":result["conservation_relative"],
                             "status":"conditional_no_slip_zero_entry_resistance"})
    assert len(rows)==144
    write_csv(OUT/"v8_水力上界逐种子.csv",rows)
    summary=[]
    l_raw=read_csv(PACKAGE/"inputs/salt/salt_transport_tidy.csv")
    l_times={}
    for r in l_raw:
        if r["source_opju"]=="Polymer_salt_2h.opju":
            continue
        if r["sample_id"].split()[0].isdigit() and int(r["sample_id"].split()[0]) in (200,350,700):
            if "dry" in r["condition"].lower() or "unwet" in r["condition"].lower():
                h=int(r["sample_id"].split()[0]); l_times[h]=max(l_times.get(h,0),float(r["time"]))
    for split in ("full_B","time_B")+tuple(f"film_{h}" for h in hs):
        group=[r for r in rows if r["split"]==split]
        for h in sorted(set(r["h_nm"] for r in group)):
            items=[r for r in group if r["h_nm"]==h]
            assert len(items)==8
            mean=np.mean([r["volume_fraction_upper_120min"] for r in items])
            summary.append({"split":split,"h_nm":h,"seed_count":8,
                            "volume_fraction_upper_120min_mean":mean,
                            "below_5pct_120min":mean<0.05,
                            "L_tmax_min_inferred":l_times.get(h,"") if split=="full_B" else "",
                            "L_volume_fraction_upper_mean":mean*l_times[h]/120 if split=="full_B" and h in l_times else "",
                            "interpretation":"conditional hydraulic bound; slip, wetting and actual pressure unknown"})
    write_csv(OUT/"v8_水力上界汇总.csv",summary)
    print(f"hydraulic solves={len(rows)}; 120min groups below 5pct={sum(r['below_5pct_120min'] for r in summary)}/{len(summary)}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=("fixed","water"))
    args=parser.parse_args()
    fixed_scenarios() if args.stage=="fixed" else water_bounds()
