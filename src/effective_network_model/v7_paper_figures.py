"""Plot v7 evidence groups 1--3 from saved tables; no solver or fitting."""

from __future__ import annotations

import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/processed/figures"
os.environ.setdefault("MPLCONFIGDIR", str(OUT))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def source(name, rows):
    with (OUT / f"v7_{name}.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["panel", "series", "h_nm", "x", "y", "metric", "source"])
        writer.writeheader()
        writer.writerows(rows)


def save(fig, name):
    fig.savefig(OUT / f"v7_{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"v7_{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def figure1():
    anchors = read("WP1/geometry_anchor_audit.csv")
    b = [r for r in read("WP4/all_predictions.csv") if r["model"] == "NS" and r["task"] == "fullB"]
    l = [r for r in read("WP4/cross_series_transfer.csv") if r["model"] == "NS"]
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    a = axes[0]
    hs = [float(r["h_nm"]) for r in anchors]
    for key, label in (("SEM_mean_nm", "fiber diameter"), ("dp_nm", "projected pore diameter")):
        ys = [float(r[key]) for r in anchors]
        a.plot(hs, ys, "o-", label=label)
        rows.extend(dict(panel="1a", series="morphology", h_nm=r["h_nm"], x=r["h_nm"], y=r[key], metric=key, source="WP1/geometry_anchor_audit.csv") for r in anchors)
    a.set(xlabel="Morphology anchor h (nm)", ylabel="Diameter (nm)", title="A. Dry morphology anchors")
    a2 = a.twinx()
    a2.plot(hs, [float(r["phi_A"]) for r in anchors], "g^:", label="projected pore area fraction")
    a2.set_ylabel("Area fraction")
    rows.extend(dict(panel="1a", series="morphology", h_nm=r["h_nm"], x=r["h_nm"], y=r["phi_A"], metric="phi_A", source="WP1/geometry_anchor_audit.csv") for r in anchors)
    handles, labels = a.get_legend_handles_labels()
    handles2, labels2 = a2.get_legend_handles_labels()
    a.legend(handles + handles2, labels + labels2, fontsize=7)
    for ax, data, panel, limit in ((axes[1], b, "1b", 120), (axes[2], l, "1c", None)):
        hs_data = sorted({int(float(r["h_nm"])) for r in data})
        for h in hs_data:
            rr = [r for r in data if int(float(r["h_nm"])) == h]
            rr.sort(key=lambda r: float(r["time_min"]))
            ax.plot([float(r["time_min"]) for r in rr], [float(r["conductivity_uS_cm"]) for r in rr], "o-", ms=2.5, lw=1, label=f"{h} nm")
            rows.extend(dict(panel=panel, series="B" if panel == "1b" else "L", h_nm=h, x=r["time_min"], y=r["conductivity_uS_cm"], metric="observed_conductivity_uS_cm", source="WP4/all_predictions.csv" if panel == "1b" else "WP4/cross_series_transfer.csv") for r in rr)
        ax.set(xlabel="Time (min)", ylabel="Conductivity (μS/cm)", title="B. Book1 six films" if panel == "1b" else "C. Long-time series")
        ax.legend(fontsize=6, ncol=2)
        if limit is not None:
            ax.set_xlim(0, limit)
    fig.suptitle("Observed data; morphology anchors are not paired with salt specimens", fontsize=11)
    fig.tight_layout()
    assert len([r for r in rows if r["panel"] == "1b"]) == 42
    assert len([r for r in rows if r["panel"] == "1c"]) == 70
    source("Fig1_material_observation", rows)
    save(fig, "Fig1_material_observation")


def figure2():
    benchmark = read("WP2/v7_straight_chain_benchmark.csv")
    rows = [dict(panel="2b", series="straight_chain_no_access", h_nm=r["h_nm"], x=r["cells"], y=r["t_lag_relative_error"], metric="lag_relative_error", source="WP2/v7_straight_chain_benchmark.csv") for r in benchmark if r["access"] == "false"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8))
    a.scatter([0, 1, 2, 3, 4], [0, 0.35, 0.65, 0.35, 0], s=[750, 340, 340, 340, 750], c=["#bdd7e7", "#74c476", "#74c476", "#74c476", "#bdd7e7"], edgecolors="black")
    ys = [0, .35, .65, .35, 0]
    for i in range(4):
        a.annotate("", xy=(i+0.84, ys[i+1]), xytext=(i+.16, ys[i]), arrowprops={"arrowstyle": "<->", "lw": 1.5})
    for x, y, label in ((0,0,"L pool"),(1,.35,"V₁K₁"),(2,.65,"V₂K₂"),(3,.35,"V₃K₃"),(4,0,"R pool")):
        a.text(x, y, label, ha="center", va="center", fontsize=8)
    a.set(xlim=(-.5,4.5), ylim=(-.35,1.15), title="A. Capacity nodes and symmetric conductances")
    a.axis("off")
    a.text(2,-.27,"ΣMᵢuᵢ = constant;  Fᵢⱼ = gᵢⱼ(uᵢ−uⱼ)",ha="center",fontsize=9)
    for h in ("200", "2600"):
        rr = [r for r in benchmark if r["access"] == "false" and r["h_nm"] == h and r["K"] == "1.0"]
        if rr:
            b.loglog([float(r["cells"]) for r in rr], [float(r["t_lag_relative_error"]) for r in rr], "o-", label=f"h={h} nm")
    b.set(xlabel="Straight-chain cells", ylabel="Continuum lag relative error", title="B. Independent straight-chain benchmark")
    b.grid(True, which="both", alpha=.2)
    b.legend(fontsize=7)
    fig.tight_layout()
    source("Fig2_model_benchmark", rows)
    save(fig, "Fig2_model_benchmark")


def figure3():
    data = read("WP2/v7_c1_main_summary.csv")
    selected = [r for r in data if r["beta"] in {"0.0", "3.0"}]
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, metric, ylabel, suffix in ((axes[0], "P_mean_m_s", "P (m/s)", "steady"), (axes[1], "u_R_120min_mean", "Right-pool u at 120 min", "finite time")):
        for series in ("H", "S"):
            for beta in ("0.0", "3.0"):
                rr = sorted([r for r in selected if r["series"] == series and r["beta"] == beta], key=lambda r: float(r["h_nm"]))
                ax.plot([float(r["h_nm"]) for r in rr], [float(r[metric]) for r in rr], "o-" if series == "H" else "s--", ms=3, label=f"{series}, β={beta}")
                rows.extend(dict(panel="3a" if suffix == "steady" else "3b", series=series, h_nm=r["h_nm"], x=r["h_nm"], y=r[metric], metric=metric, source="WP2/v7_c1_main_summary.csv") for r in rr)
        ax.set(xscale="log", yscale="log", xlabel="h (nm)", ylabel=ylabel, title=f"{suffix.capitalize()} comparison; eight-seed means")
        ax.legend(fontsize=7)
        ax.grid(True, which="both", alpha=.2)
    fig.suptitle("H fixes draw357 morphology; S varies morphology conditionally", fontsize=10)
    fig.tight_layout()
    source("Fig3_structure_vs_thickness", rows)
    save(fig, "Fig3_structure_vs_thickness")


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    print("v7 paper figures 1-3: saved SVG, PNG and source CSV; no model rerun")
