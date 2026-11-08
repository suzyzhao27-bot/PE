"""Render review figures from fixed A and D tables; no model recomputation."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
A = DATA / "25_论文裁剪A厚薄保守界_20261031.csv"
D = DATA / "25_论文裁剪D尿素浓度条件范围_20261031.csv"
D_GAIN = DATA / "25_论文裁剪D尿素处理响应范围_20261031.csv"
D_GAIN_SENS = DATA / "25_论文裁剪D迁移系数敏感性_处理响应_20261031.csv"
OUT_A = DATA / "25_论文裁剪A主图预览_20261031.png"
OUT_D = DATA / "25_论文裁剪D主图预览_20261031.png"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def a_figure() -> None:
    rows = read(A)
    groups = (
        ("G0 E1\nδ=0", "G0", 0, "E1", 1, 0),
        ("G0 E1\nδ=.20", "G0", .2, "E1", 1, 0),
        ("G1 E1\nδ=.20", "G1", .2, "E1", 1, 0),
        ("G1 E2\nΓ=10", "G1", .2, "E2", 10, 0),
        ("G1 E1\n±5% Dbar", "G1", .2, "E1", 1, .05),
        ("G2 E1\nδ=.20", "G2", .2, "E1", 1, 0),
    )
    fig, ax = plt.subplots(figsize=(10.7, 4.7))
    colors = {"C1_T85": "#1769aa", "C2_T85": "#c65b25"}
    for scenario, offset in (("C1_T85", -.11), ("C2_T85", .11)):
        values = []
        for _, family, delta, external, gamma, sensitivity in groups:
            selected = [r for r in rows if r["scenario"] == scenario and r["window"] == "W2"
                        and r["family"] == family and float(r["delta_rho"]) == delta
                        and r["external"] == external and int(r["Gamma"]) == gamma
                        and float(r["mobility_sensitivity_fraction"]) == sensitivity]
            if len(selected) != 1:
                raise ValueError("A figure condition identity mismatch")
            values.append(float(selected[0]["lower_b_B2600_over_B65"]))
        ax.scatter(np.arange(len(groups)) + offset, values, color=colors[scenario],
                   s=65, label=scenario[:2], zorder=3)
    ax.axhline(1, color="#555", lw=1, ls="--")
    ax.set_yscale("log")
    ax.set_ylim(.8, 6000)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups])
    ax.set_ylabel("Necessary lower bound on b(B2600) / b(B65)")
    ax.set_title("Archived-model conditional bounds, W2 (60–100 min)")
    ax.legend(title="Conductivity conversion", frameon=False)
    ax.grid(axis="y", which="both", alpha=.22)
    fig.text(.015, .048, "Except G0 δ=0, all groups use δ=.20; G0/G1 also differ in geometry and charge range. Γ=10 and ±5% are scenarios.",
             fontsize=8, color="#555")
    fig.text(.015, .012, "Cancellation: independent external resistance per curve removes data-derived positive bounds (one window).",
             fontsize=8, color="#555")
    fig.tight_layout(rect=(0, .085, 1, 1))
    fig.savefig(OUT_A, dpi=190)
    plt.close(fig)


def d_figure() -> None:
    rows = [r for r in read(D) if r["external_relation"] == "DIFFUSION_SCALED_SHARED"]
    samples = ("B65", "R65", "B700", "R700", "B2600", "R2600")
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.0),
                             gridspec_kw={"width_ratios": (1.15, 1.15, 1)})
    axes[1].sharey(axes[0])
    colors = {"C1_T85": "#1769aa", "C2_T85": "#c65b25"}
    for ax, time in zip(axes[:2], (60, 100)):
        for scenario, offset in (("C1_T85", -.12), ("C2_T85", .12)):
            selected = []
            for sample in samples:
                found = [r for r in rows if r["scenario"] == scenario
                         and r["sample_id"] == sample and int(r["time_min"]) == time]
                if len(found) != 1:
                    raise ValueError("D figure condition identity mismatch")
                selected.append(found[0])
            lo = np.array([float(r["urea_receiver_min_mM"]) for r in selected])
            hi = np.array([float(r["urea_receiver_max_mM"]) for r in selected])
            x = np.arange(len(samples)) + offset
            ax.vlines(x, lo, hi, lw=2, color=colors[scenario], label=scenario[:2])
            ax.hlines(lo, x - .045, x + .045, lw=1.5, color=colors[scenario])
            ax.hlines(hi, x - .045, x + .045, lw=1.5, color=colors[scenario])
        ax.set_xticks(range(len(samples)), samples)
        ax.set_yscale("log")
        ax.set_title(f"Old-device conditional target, {time} min")
        ax.grid(axis="y", which="both", alpha=.22)
    axes[0].set_ylabel("Predicted urea receiver concentration (mM)")
    handles, labels = axes[1].get_legend_handles_labels()
    axes[1].legend(handles[:2], labels[:2], title="Conversion", frameon=False,
                   loc="lower right")
    base = [r for r in read(D_GAIN) if r["external_relation"] == "DIFFUSION_SCALED_SHARED"]
    sens = read(D_GAIN_SENS)
    ax = axes[2]
    for scenario, offset in (("C1_T85", -.11), ("C2_T85", .11)):
        for x0, h in enumerate((65, 700, 2600)):
            old = [r for r in base if r["scenario"] == scenario and r["h_nm"] == str(h)
                   and r["window"] == "W2"]
            wide = [r for r in sens if r["scenario"] == scenario and r["h_nm"] == str(h)
                    and r["window"] == "W2"]
            if len(old) != 1 or len(wide) != 1:
                raise ValueError("D gain figure identity mismatch")
            x = x0 + offset
            ax.vlines(x, float(wide[0]["mobility_5pct_min"]),
                      float(wide[0]["mobility_5pct_max"]),
                      lw=7, alpha=.25, color=colors[scenario])
            ax.vlines(x, float(old[0]["urea_gain_min"]),
                      float(old[0]["urea_gain_max"]),
                      lw=2, color=colors[scenario])
            ax.scatter(x, float(old[0]["salt_gain_R_over_B"]), marker="D", s=21,
                       facecolor="white", edgecolor="#222", linewidth=.8, zorder=4)
    ax.axhline(1, color="#444", lw=1, ls="--")
    ax.set_yscale("log")
    ax.set_ylim(.65, 500)
    ax.set_xticks((0, 1, 2), ("65", "700", "2600"))
    ax.set_xlabel("Nominal thickness (nm)")
    ax.set_ylabel("Predicted urea R/B rate ratio")
    ax.set_title("W2 model response; unity reference")
    ax.grid(axis="y", which="both", alpha=.22)
    fig.text(.015, .048, "Left: equal-volume, no-sampling/no-net-flow model. Right: legacy salt ratios and urea main / ±5% ranges.",
             fontsize=8, color="#555")
    fig.text(.015, .012, "Cancellation: independent urea external resistance removes finite response bounds; no cancellation value is plotted.",
             fontsize=8, color="#555")
    fig.tight_layout(rect=(0, .085, 1, 1))
    fig.savefig(OUT_D, dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    a_figure()
    d_figure()
    print(f"figures={OUT_A.name},{OUT_D.name}")
