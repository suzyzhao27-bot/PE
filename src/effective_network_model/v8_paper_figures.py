"""Draw v8.0.1 evidence panels from saved tables; never fit while plotting."""

import csv
import os
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR",str(ROOT/".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT=(ROOT/"figures") if (ROOT/"taskbook").is_dir() else ROOT/"data/processed/figures"


def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def save(fig,name):
    fig.tight_layout()
    fig.savefig(OUT/f"v8_{name}.png",dpi=180)
    fig.savefig(OUT/f"v8_{name}.svg")
    plt.close(fig)


def fig1():
    source=rows(ROOT/"WP4/v8_全部逐点预测.csv")
    observed=[r for r in source if r["model"]=="Tcal" and r["split"]=="full_B"]
    fig,ax=plt.subplots(figsize=(7.4,4.6))
    for h in (65,120,200,350,700,2600):
        data=sorted([r for r in observed if int(r["h_nm"])==h],key=lambda r:float(r["t_min"]))
        ax.plot([float(r["t_min"]) for r in data],[float(r["delta_obs_uScm"]) for r in data],
                marker="o",markersize=3,label=f"{h} nm")
    ax.set_yscale("symlog",linthresh=1)
    ax.set(xlabel="Time (min; Book1 explicit units)",ylabel="Added conductivity (uS/cm)",
           title="B series: original added conductivity")
    ax.legend(ncol=3,fontsize=8)
    save(fig,"Fig1_B_observations")


def fig2():
    data=rows(ROOT/"WP1/量纲与假设表.csv")
    h=np.array([int(r["h_nm"]) for r in data])
    lag=np.array([float(r["lag_upper_alpha10_epsilon035_s"]) for r in data])
    chi=np.array([float(r["chi_app_alpha10"]) for r in data])
    fig,(a,b)=plt.subplots(1,2,figsize=(9.4,3.9))
    a.bar(range(len(h)),lag,color="#446f91")
    a.axhline(1200,color="#b33939",ls="--",label="20 min")
    a.set_yscale("log");a.set_xticks(range(len(h)),h)
    a.set(xlabel="Thickness (nm)",ylabel="Lag upper bound (s)",title="Conditional: alpha=10, K<=1")
    a.legend(fontsize=8)
    b.plot(h,chi/chi[3],marker="o",color="#b76631")
    b.set_yscale("log");b.set_xscale("log")
    b.set(xlabel="Thickness (nm)",ylabel="Apparent chi / chi(350 nm)",title="Equal-area endpoint requirement")
    save(fig,"Fig2_constraints")


def fig3():
    data=rows(ROOT/"WP2/预筛全部记录.csv")
    sample=sorted([r for r in data if int(r["h_nm"])==2600 and int(r["seed"])==7001],
                  key=lambda r:float(r["q"]))
    q=np.array([float(r["q"]) for r in sample])
    upper=np.array([float(r["seed_group_mean_uScm"]) for r in sample])
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    ax.plot(q,upper,marker="o",ms=3,color="#365d8d")
    ax.axhline(1363.25/3,color="#bb3838",ls="--",label="2600 nm observed / 3")
    ax.axvline(0.02813067042841564,color="#555555",ls=":",label="first not excluded q")
    ax.set_xscale("log");ax.set_yscale("symlog",linthresh=1)
    ax.set_ylim(bottom=0)
    ax.set(xlabel="Common throat factor q",ylabel="Mean upper endpoint (uS/cm)",
           title="Necessary bound, not a network prediction")
    ax.legend(fontsize=8)
    save(fig,"Fig3_prescreen")


def fig4():
    data=rows(ROOT/"WP4/v8_全部逐点预测.csv")
    full=[r for r in data if r["split"]=="full_B"]
    fig,axes=plt.subplots(2,3,figsize=(11.4,6.4),sharex=True)
    for ax,h in zip(axes.flat,(65,120,200,350,700,2600)):
        reference=sorted([r for r in full if r["model"]=="Tcal" and int(r["h_nm"])==h],
                         key=lambda r:float(r["t_min"]))
        ax.scatter([float(r["t_min"]) for r in reference],
                   [float(r["delta_obs_uScm"]) for r in reference],color="black",s=14,label="observed")
        for model,color in (("Ncal","#9d9d9d"),("Scal","#3a77aa"),("Tcal","#c4643b")):
            part=sorted([r for r in full if r["model"]==model and int(r["h_nm"])==h],
                        key=lambda r:float(r["t_min"]))
            ax.plot([float(r["t_min"]) for r in part],[float(r["delta_pred_uScm"]) for r in part],
                    color=color,label=model,lw=1.3)
        ax.set_yscale("symlog",linthresh=1)
        ax.set_title(f"{h} nm",fontsize=10)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Added conductivity (uS/cm)")
    axes.flat[0].legend(fontsize=7,ncol=2)
    save(fig,"Fig4_fullB_comparison")


def fig5():
    data=rows(ROOT/"WP4/v8_模型比较.csv")
    models=[r["model"] for r in data]
    x=np.arange(3)
    fig,(a,b)=plt.subplots(1,2,figsize=(9.2,3.8))
    a.bar(x-0.16,[float(r["E_time"]) for r in data],width=.32,label="time holdout")
    a.bar(x+0.16,[float(r["E_film"]) for r in data],width=.32,label="film holdout")
    a.axhline(.7,color="#ba3434",ls="--",label="frozen 0.70 gate")
    a.set_xticks(x,models);a.set(ylabel="Equal-film asinh E",title="Retrospective holdouts")
    a.legend(fontsize=7)
    b.bar(x,[float(r["film_2600_endpoint_pred_over_obs"]) for r in data],color="#a76a34")
    b.axhspan(1/3,3,color="#559966",alpha=.2,label="frozen 1/3 to 3 gate")
    b.set_yscale("log");b.set_xticks(x,models)
    b.set(ylabel="2600 nm predicted / observed",title="Whole-film endpoint")
    b.legend(fontsize=7)
    save(fig,"Fig5_prediction")


def fig6():
    fixed=rows(ROOT/"WP5/v8_固定十情景汇总.csv")
    domain=rows(ROOT/"WP5/v8_域尺寸汇总.csv")
    fig,(a,b)=plt.subplots(1,2,figsize=(10,4.4))
    labels=[r["scenario"] for r in fixed]
    a.barh(labels,[float(r["film_2600_endpoint_pred_over_obs"]) for r in fixed],color="#9c7044")
    a.axvspan(1/3,3,color="#559966",alpha=.2)
    a.set_xscale("log");a.set(xlabel="2600 nm predicted / observed",title="Ten fixed scenarios")
    d=sorted(domain,key=lambda r:float(r["relative_mean_change"]),reverse=True)
    labels=[f"{r['model']} {r['split']} {r['h_nm']}" for r in d[:8]]
    b.barh(labels,[100*float(r["relative_mean_change"]) for r in d[:8]],color="#487d9e")
    b.axvline(10,color="#b33939",ls="--")
    b.set(xlabel="12 to 24 mean P change (%)",title="Largest domain effects")
    save(fig,"Fig6_robustness")


if __name__=="__main__":
    for fn in (fig1,fig2,fig3,fig4,fig5,fig6):
        fn()
    print("v8 figures: 6 SVG + 6 PNG from saved tables")
