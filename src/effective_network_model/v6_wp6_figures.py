"""Build three summary figures from frozen WP4/WP5 tables, without refitting."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP6" / "figures"
os.environ.setdefault("MPLCONFIGDIR", str(OUT))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    OUT.mkdir(exist_ok=True)
    colors = {"M0": "#567a93", "M1": "#d17a3b", "M2": "#6f9b66"}
    residuals = read("WP4/03_残差与模型比较.csv")
    hs = [65, 120, 200, 350, 700, 2600]
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), layout="constrained")
    for panel, key, title in zip(ax, ("rmse_early_asinh", "rmse_late_asinh"),
                                 ("Training: <=60 min", "Training: >60 min")):
        for j, model in enumerate(("M0", "M1", "M2")):
            values = [float(next(r[key] for r in residuals if r["model_id"] == model and
                                 int(r["thickness_nm"]) == h)) for h in hs]
            panel.bar(np.arange(6) + (j-1)*0.24, values, width=0.23,
                      label=model, color=colors[model])
        panel.set_xticks(np.arange(6), [str(h) for h in hs])
        panel.set_xlabel("Thickness (nm)")
        panel.set_ylabel("RMSE in asinh units")
        panel.set_title(title)
        panel.grid(axis="y", alpha=0.2)
    ax[0].legend(frameon=False, ncol=3)
    fig.savefig(OUT / "04_训练残差的逐膜分解.png", dpi=180)
    plt.close(fig)

    profiles = read("WP5/03_参数辨识与稳健性结果.csv")
    panels = [("M1", "log10_gamma_per_m"), ("M1", "log10_q_conductivity"),
              ("M2", "log10_chi"), ("M2", "log10_tau_s")]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.4), layout="constrained")
    for axis, (model, parameter) in zip(axes.flat, panels):
        rows = [r for r in profiles if r["model_id"] == model and
                r["fixed_parameter"] == parameter]
        xs = [float(r["fixed_value_transformed"]) for r in rows]
        ys = [float(r["loss"]) if r["loss"] and r["success"] == "True" else np.nan
              for r in rows]
        axis.plot(xs, ys, "o-", color=colors[model], markersize=3)
        if any(r["success"] != "True" for r in rows):
            axis.plot([float(r["fixed_value_transformed"]) for r in rows
                       if r["success"] != "True"],
                      [min(x for x in ys if np.isfinite(x))]*sum(r["success"] != "True" for r in rows),
                      "x", color="crimson", label="failed")
            axis.legend(frameon=False)
        axis.set_title(f"{model}: {parameter}")
        axis.set_xlabel("Fixed log10 parameter")
        axis.set_ylabel("Reoptimized training loss")
        axis.grid(alpha=0.2)
    fig.savefig(OUT / "05_共享参数条件代价剖面.png", dpi=180)
    plt.close(fig)

    sensitivity = read("WP5/07_固定情景与读数扰动.csv")
    base = float(json.loads((ROOT / "WP4" / "M1_拟合摘要.json").read_text(encoding="utf-8"))["best_loss"])
    fixed = [r for r in sensitivity if not r["scenario_id"].startswith("reading_")]
    fixed.sort(key=lambda r: abs(float(r["loss_bdf"])/base-1), reverse=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
    for ax, subset, title in ((axes[0], fixed[:2], "Largest changes"),
                              (axes[1], fixed[2:], "Remaining fixed scenarios")):
        x = [r["scenario_id"] for r in subset]
        y = [float(r["loss_bdf"])/base-1 for r in subset]
        ax.barh(x[::-1], y[::-1], color="#8e77a9")
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Loss change / frozen M1")
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.2)
    fig.savefig(OUT / "06_固定情景的训练损失变化.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
