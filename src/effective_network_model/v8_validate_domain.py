"""Frozen v8.0.1 section 5.2: selected-point 12-to-24 box comparison."""

import math
from pathlib import Path

import numpy as np

from v7_network_geometry import generate, morphology
from v8_compare import read_csv, setup, write_csv
from v8_transport import steady, weighted_graph

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"WP5"


def main():
    new,old=setup()
    hs=new["geometry"]["B_h_nm"]
    seeds=new["geometry"]["main_seeds"]
    selected={(r["model"],r["split"]):r for r in read_csv(ROOT/"WP4/v8_最优参数.csv")}
    base=read_csv(ROOT/"WP4/v8_稳态表.csv")
    base_map={(r["model"],int(r["material_index"]),int(r["h_nm"]),int(r["seed"])):float(r["P_m_s"])
              for r in base}
    old["geometry"].update(nx=24,ny=24)
    max_nz=max(math.ceil(h/sum(morphology(h,"S",old["geometry"]["anchors"])[:2])) for h in hs)
    rows=[]
    for h in hs:
        tasks=[(model,"full_B") for model in ("Ncal","Scal","Tcal")]
        tasks.extend(("Tcal",split) for split in ("time_B",f"film_{h}"))
        for seed in seeds:
            graph,meta=generate(h,"S",seed,old,max_nz)
            info={"s_nm":meta["s_nm"],"nz":meta["nz"],"epsilon":meta["epsilon"],"nx":24,"ny":24}
            for model,split in tasks:
                choice=selected[(model,split)]
                q,beta=float(choice["q"]),float(choice["beta"])
                local="Tcal" if model=="Tcal" else "Scal"
                edges,area=weighted_graph(graph,info,q,local,beta)
                solved=steady(len(graph["nodes"]),edges,area,rtol=1e-9,solver="cg")
                key=(local,0 if model=="Ncal" else int(choice["material_index"]),h,seed)
                # Tcal indices in selected table enumerate only q values not excluded by G2.
                base_p=base_map[key]
                rows.append({"model":model,"split":split,"h_nm":h,"seed":seed,
                             "q":q,"beta":beta,"P_12_m_s":base_p,
                             "P_24_m_s":solved["P_m_s"],
                             "relative_change_vs_12":abs(solved["P_m_s"]-base_p)/base_p if base_p else "",
                             "conservation_relative_24":solved["conservation_relative"],
                             "solver":"CG Jacobi rtol=1e-9","status":"PASS_TECHNICAL"})
    assert len(rows)==240
    write_csv(OUT/"v8_域尺寸逐种子.csv",rows)
    grouped={}
    for row in rows:
        grouped.setdefault((row["model"],row["split"],row["h_nm"]),[]).append(row)
    summary=[]
    for (model,split,h),items in grouped.items():
        assert len(items)==8
        p12=np.mean([r["P_12_m_s"] for r in items])
        p24=np.mean([r["P_24_m_s"] for r in items])
        difference=abs(p24/p12-1) if p12 else ""
        summary.append({"model":model,"split":split,"h_nm":h,"seed_count":8,
                        "P_12_mean_m_s":p12,"P_24_mean_m_s":p24,
                        "relative_mean_change":difference,
                        "within_10pct":difference<=0.1 if p12 else False,
                        "interpretation":"same seed labels, different 12/24 random-array shapes; compare distributions"})
    write_csv(OUT/"v8_域尺寸汇总.csv",summary)
    print(f"domain solves={len(rows)}; 10pct pass={sum(r['within_10pct'] for r in summary)}/{len(summary)}")


if __name__=="__main__":
    main()
