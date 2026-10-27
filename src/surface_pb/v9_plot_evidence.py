"""Plot v9 evidence from saved CSV only; no field solve or refit."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "data/processed/figures"


def read_csv(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def main():
    calibration = read_csv("WP4/v9_校准七档汇总.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    gamma = np.array([float(r["gamma_nm2"]) for r in calibration])
    for index, color in ((1, "#1f77b4"), (2, "#d62728")):
        y = np.array([float(r[f"start{index}_E_cal"]) for r in calibration])
        ax.plot(gamma, y, "-", color=color, lw=1, alpha=0.65)
        for x, value, row in zip(gamma, y, calibration):
            stable = row[f"start{index}_status"] == "STABLE_LOCAL"
            ax.plot(x, value, marker="o" if stable else "s", ms=6,
                    color=color, markerfacecolor=color if stable else "white")
    ax.axhline(1.0, ls="--", color="black", lw=1, label="Frozen compatibility gate E=1")
    ax.annotate("minimum 1.0037188 (fails gate)", (0.05, 1.0037188),
                xytext=(0.057, 2.0), arrowprops={"arrowstyle": "-", "lw": 0.8}, fontsize=8)
    ax.set_xscale("log")
    ax.set_xticks(gamma, [f"{x:g}" for x in gamma])
    ax.set_ylim(0, 11)
    ax.set_xlabel(r"Fixed site capacity $\Gamma$ (nm$^{-2}$)")
    ax.set_ylabel(r"Training $E_{cal}$ on 182 pH 3/7/10")
    ax.set_title("Frozen calibration; no compatible parameter set")
    ax.grid(axis="y", alpha=0.25)
    ax.text(0.98, 0.96, "MODEL NOT MATERIAL-VALIDATED\nFilled circle: stable local; open square: unresolved",
            ha="right", va="top", transform=ax.transAxes, fontsize=7)
    fig.tight_layout()
    save(fig, "v9_calibration_gate")

    mapped = read_csv("WP5/v9_MODEL_ONLY_精简厚度图.csv")
    errors = read_csv("WP5/v9_MODEL_ONLY_N4_N11数值误差.csv")
    sensitivity = read_csv("WP5/v9_MODEL_ONLY_八项敏感性差值.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0))
    ax = axes[0]
    colors = {0.01: "#9467bd", 0.1: "#1f77b4", 1.0: "#2ca02c", 10.0: "#ff7f0e"}
    for salt, color in colors.items():
        relevant = sorted((r for r in mapped if r["kind"] == "free"
                           and float(r["salt_mM"]) == salt), key=lambda r: int(r["N"]))
        ax.plot([int(r["N"]) for r in relevant],
                [abs(float(r["sigma_mean_C_m2"])) for r in relevant],
                marker="o", lw=1.3, color=color, label=f"{salt:g} mM")
    by_N = {N: max(float(r["absolute_change_C_m2"]) for r in errors if int(r["N"]) == N)
            for N in (4, 11)}
    sigma1 = {int(r["N"]): abs(float(r["sigma_mean_C_m2"])) for r in mapped
              if r["kind"] == "free" and float(r["salt_mM"]) == 1}
    ax.errorbar([4, 11], [sigma1[4], sigma1[11]], yerr=[by_N[4], by_N[11]],
                fmt="none", ecolor="black", capsize=3, lw=1, label="1 mM mesh/domain change")
    ax.set_xlabel("Number of fiber layers N")
    ax.set_ylabel(r"Free-membrane $|\sigma_{mean}|$ (C m$^{-2}$)")
    ax.set_ylim(0, 0.009)
    ax.set_title("Pure thickness, fixed d/g, pH 7")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.2)
    ax = axes[1]
    order = ["g_x0.5", "g_x2", "H_x0.9", "H_x1.1", "cs_x0.8", "cs_x1.2",
             "PE_eps_2.3", "aligned_layers"]
    values = [max(abs(float(r["delta_zeta_mV"])) for r in sensitivity
                  if r["variation"] == variation and r["kind"] == "backed") for variation in order]
    x = np.arange(len(order))
    ax.bar(x, values, color=["#2c7fb8" if k != "PE_eps_2.3" else "#7b3294" for k in order])
    ax.axhline(5, ls="--", color="black", lw=1, label="Frozen dielectric risk gate: 5 mV")
    ax.set_xticks(x, ["g/2", "2g", "0.9H", "1.1H", "0.8c", "1.2c", "finite PE", "aligned"],
                  rotation=45, ha="right")
    ax.set_ylabel(r"Maximum $|\Delta\zeta_{app}|$ across nine conditions (mV)")
    ax.set_ylim(0, 5.5)
    ax.set_title("One-factor checks; no refit")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.suptitle("MODEL ONLY — legacy chemistry; not PE material validation", fontsize=10)
    fig.tight_layout()
    save(fig, "v9_model_only_thickness_sensitivity")
    print("v9 figures: 2 SVG + 2 PNG from saved CSV")


if __name__ == "__main__":
    main()
