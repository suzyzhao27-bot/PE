"""Render first-round manuscript panel layouts from fixed design and archived A data.

No experimental points are synthesized. The only numerical panel reads the archived
conditional A table. Existing A/D review plots and original source files stay intact.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "salt_transport"
A_TABLE = OUT / "25_论文裁剪A厚薄保守界_20261031.csv"
AUDIT = Path(r"C:\Users\12989\.agents\skills\nature-figure\scripts")
sys.path.insert(0, str(AUDIT))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

INK = "#22364a"
BLUE = "#2673a6"
ORANGE = "#bd6334"
MUTED = "#627181"
PALE = "#f1f4f6"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": "white",
    }
)


def grid(number: int, nrows: int, ncols: int, titles: list[str]):
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.2, 2.2 * nrows + 0.65))
    flat = np.asarray(axes).ravel()
    for idx, ax in enumerate(flat):
        if idx >= len(titles):
            ax.set_visible(False)
            continue
        ax.set_title(f"{chr(97 + idx)}  {titles[idx]}", loc="left", fontsize=7.5,
                     fontweight="bold", color=INK, pad=7)
    fig.suptitle(f"Figure {number} | First-round panel architecture", x=0.05,
                 ha="left", fontsize=9, fontweight="bold", color=INK)
    return fig, flat


def pending(ax, lines: str):
    ax.set_facecolor(PALE)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("#b9c4cc")
        spine.set_linestyle("--")
    ax.text(0.5, 0.58, "DATA PENDING", ha="center", va="center",
            transform=ax.transAxes, color=MUTED, fontsize=8, fontweight="bold")
    ax.text(0.5, 0.37, lines, ha="center", va="center", transform=ax.transAxes,
            color=MUTED, fontsize=8, linespacing=1.5)


def design(ax, text: str):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(0.5, 0.5, text, ha="center", va="center", color=INK,
            fontsize=8, linespacing=1.6,
            bbox={"facecolor": "#e9f1f6", "edgecolor": BLUE, "boxstyle": "round,pad=0.7"})


def save(fig, number: int):
    fig.subplots_adjust(left=0.16, right=0.97, top=0.84, bottom=0.14,
                        hspace=0.62, wspace=0.52)
    fig.canvas.draw()
    prefix = OUT / f"20261102_first_round_Fig{number}"
    require_matplotlib_panel_alignment(
        fig, json_out=str(prefix) + ".alignment.json", tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, strict=True
    )
    fig.savefig(str(prefix) + ".pdf")
    fig.savefig(str(prefix) + ".svg")
    fig.savefig(str(prefix) + ".png", dpi=360)
    plt.close(fig)


def fig1():
    fig, ax = grid(1, 3, 2, ["Fiber network", "Paired H-cell design",
                                "Measured thickness", "High-salt NaCl histories",
                                "Parent-paired salt gain"])
    pending(ax[0], "Source-verified image\nand calibrated scale bar")
    design(ax[1], "One parent membrane → B / E pieces\nEach run: NaCl + urea\nRecord volumes and samples")
    pending(ax[2], "Nine parent membranes;\nmeasured h and identity")
    pending(ax[3], "18 high-salt runs; 0–240 min\nconcentration + sampling events")
    pending(ax[4], "Three parent pairs per thickness\nE/B salt gain + interval")
    save(fig, 1)


def fig2():
    fig, ax = grid(2, 2, 2, ["External resistance", "Urea perturbation screen",
                                "Water waiting versus prewetting", "One versus three layers"])
    pending(ax[0], "Six no-PE controls;\nboth probes and matched membrane runs")
    pending(ax[1], "Four P1 screen runs;\nwith / without urea")
    pending(ax[2], "B / W / E matched conditions;\nwater waiting control")
    pending(ax[3], "One / three E layers;\ntransfer + resistance")
    save(fig, 2)


def fig3():
    fig, ax = grid(3, 1, 3, ["Within-run dual probe", "Paired common response",
                                "Ion-specific remainder"])
    pending(ax[0], "One parent membrane;\nNaCl + urea per run")
    pending(ax[1], "log E_s vs log E_u;\nthree pairs / thickness")
    pending(ax[2], "Π = E_s / E_u;\nlow vs high salt")
    save(fig, 3)


def archived_a(ax):
    with A_TABLE.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    groups = [
        ("G0, δ=0", "G0", 0.0, "E1", 1, 0.0),
        ("G0, δ=.20", "G0", .2, "E1", 1, 0.0),
        ("G1, δ=.20", "G1", .2, "E1", 1, 0.0),
        ("G1, Γ=10", "G1", .2, "E2", 10, 0.0),
        ("G1, ±5%", "G1", .2, "E1", 1, .05),
        ("G2, δ=.20", "G2", .2, "E1", 1, 0.0),
    ]
    for scenario, color, shift in [("C1_T85", BLUE, -.11), ("C2_T85", ORANGE, .11)]:
        values = []
        for _, family, delta, external, gamma, sensitivity in groups:
            hits = [r for r in rows if r["scenario"] == scenario and r["window"] == "W2"
                    and r["family"] == family and float(r["delta_rho"]) == delta
                    and r["external"] == external and int(r["Gamma"]) == gamma
                    and float(r["mobility_sensitivity_fraction"]) == sensitivity]
            if len(hits) != 1:
                raise ValueError("Archived A scenario identity mismatch")
            values.append(float(hits[0]["lower_b_B2600_over_B65"]))
        ax.scatter(values, np.arange(len(groups)) + shift, color=color, s=13,
                   label=scenario[:2], zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(range(len(groups)), [g[0] for g in groups])
    ax.invert_yaxis()
    ax.axvline(1, color=MUTED, lw=.8, ls="--")
    ax.set_xlabel("Necessary lower bound, b(B2600)/b(B65)")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper left")
    ax.grid(axis="x", alpha=.2)


def fig4():
    fig, ax = grid(4, 2, 2, ["Series-resistance model", "Urea-constrained capacity",
                                "Path-length-only reference", "Archived necessary bounds"])
    design(ax[0], "Device resistance + pore-mouth resistance\n+ membrane path resistance (h, a, b)\nb = εf/τ: composite capacity, not pore count")
    pending(ax[1], "Measured urea + external resistance\n→ independent b interval")
    h = np.array([65, 700, 2600], dtype=float)
    ax[2].plot(h, h / h[0], color=BLUE, lw=1.7, marker="o", ms=3)
    ax[2].set_xscale("log")
    ax[2].set_yscale("log")
    ax[2].set_xticks(h, ["65", "700", "2600"])
    ax[2].set_xlabel("Nominal path length h (nm)")
    ax[2].set_ylabel("Path resistance / 65 nm reference")
    ax[2].grid(alpha=.2)
    archived_a(ax[3])
    fig.text(.16, .035, "4c: fixed-capacity path term only. 4d: legacy W2, shared external resistance; T12 unresolved; no new-batch data.",
             fontsize=8, color=MUTED)
    save(fig, 4)


def fig5():
    fig, ax = grid(5, 2, 2, ["Calibration and held-out design",
                                "Predicted versus observed salt",
                                "All held-out residuals", "Discrete evidence status"])
    matrix = [
        ["CAL", "CAL", "HOLD", "HOLD"],
        ["—", "—", "HOLD", "HOLD"],
        ["HOLD", "HOLD", "HOLD", "HOLD"],
    ]
    colors = {"CAL": "#d7e8f3", "HOLD": "#f6e7da", "—": PALE}
    for row in range(3):
        for col in range(4):
            value = matrix[row][col]
            ax[0].add_patch(plt.Rectangle((col, row), .98, .98, facecolor=colors[value],
                                           edgecolor="white"))
            ax[0].text(col + .49, row + .49, value, ha="center", va="center",
                       color=INK, fontsize=8, fontweight="bold")
    ax[0].set_xlim(0, 4)
    ax[0].set_ylim(3, 0)
    ax[0].set_xticks(np.arange(4) + .5, ["Low\nB", "Low\nE", "High\nB", "High\nE"])
    ax[0].set_yticks(np.arange(3) + .5, ["Thin", "Middle", "Thick"])
    ax[0].tick_params(length=0)
    pending(ax[1], "Eight held-out conditions;\nurea-constrained b; measured salt")
    pending(ax[2], "All eight signed residuals;\nneutral + charged references")
    pending(ax[3], "3 thicknesses × 2 salinities;\ncommon / mixed / ion-extra / unresolved")
    save(fig, 5)


if __name__ == "__main__":
    for number, func in enumerate((fig1, fig2, fig3, fig4, fig5), 1):
        func()
        print(f"Figure {number}: PDF and PNG saved")
