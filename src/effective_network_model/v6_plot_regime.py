"""Plot the fixed WP3 synthetic finite-time slope matrix from its CSV."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "WP3/02_消融及适用域结果.csv"
TARGET = ROOT / "WP3/figures/有限时间厚度斜率适用图.png"


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r["record_type"] == "regime_map"]
    scenarios = [("M1", 0., 1.), ("M2", 10., 1.), ("M2", 1000., 1.),
                 ("M1", 0., 1e-8), ("M2", 10., 1e-8), ("M2", 1000., 1e-8)]
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.2), sharex=True, sharey=True,
                             constrained_layout=True)
    ell = 900e-9
    db = 1.6070833333333333e-9
    for ax, (model, chi, ratio) in zip(axes.flat, scenarios):
        selected = [r for r in rows if r["model_id"] == model
                    and abs(float(r["chi"]) - chi) < 1e-9
                    and abs(float(r["dm_over_db"]) / ratio - 1) < 1e-10]
        assert len(selected) == 41 * 41
        etas = sorted({float(r["eta"]) for r in selected})
        omegas = sorted({float(r["omega"]) for r in selected})
        lookup = {(float(r["eta"]), float(r["omega"])): r for r in selected}
        slopes = np.array([[float(lookup[(e, o)]["S_T"]) if lookup[(e, o)]["S_T"] else np.nan
                            for e in etas] for o in omegas])
        mesh = ax.pcolormesh(etas, omegas, np.ma.masked_invalid(np.clip(slopes, -3, 3)),
                             shading="nearest", cmap="RdBu_r", vmin=-3, vmax=3)
        clean = np.ma.masked_invalid(slopes)
        ax.contour(etas, omegas, clean, levels=[0], colors="black", linewidths=.8)
        omega_120 = ratio * db * 7200 / (ell * ell)
        if omegas[0] <= omega_120 <= omegas[-1]:
            ax.axhline(omega_120, color="gray", ls="--", lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{model}, $D_m/D_b={ratio:g}$" + (f", $\\chi={chi:g}$" if model == "M2" else ""))
        ax.grid(alpha=.15, which="both")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\eta=h/\ell$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\Omega=D_mT/\ell^2$")
    fig.colorbar(mesh, ax=axes, shrink=.82, label=r"$S_T=\partial_{\ln h}\ln M_R(T)$ (clipped to [-3,3])")
    fig.suptitle("Synthetic fixed-parameter finite-time thickness slopes; white = unresolved response",
                 fontsize=13)
    fig.savefig(TARGET, dpi=220)
    plt.close(fig)
    print(TARGET)


if __name__ == "__main__":
    main()
