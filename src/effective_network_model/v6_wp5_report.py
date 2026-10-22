"""Merge frozen WP5 holdouts and draw transparent comparison figures."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP5"
os.environ.setdefault("MPLCONFIGDIR", str(OUT))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read(name):
    with (OUT / name).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    models = ("M0", "M1", "M2")
    thickness = sum((read(f"{m}_厚度留出.csv") for m in models), [])
    time = sum((read(f"{m}_早段预测.csv") for m in models), [])
    write("01_厚度留出结果.csv", thickness)
    write("02_早段预测后段结果.csv", time)
    hvalues = (65, 120, 200, 350, 700, 2600)
    metrics = {}
    for m in models:
        h = [r for r in thickness if r["model_id"] == m]
        t = [r for r in time if r["model_id"] == m]
        metrics[m] = {"E_h_macro": float(np.mean([float(r["E_h"]) for r in h])),
                      "E_time_macro": float(np.mean([float(r["E_time"]) for r in t])),
                      "per_film_h": {r["thickness_nm"]: float(r["E_h"]) for r in h},
                      "per_film_time": {r["thickness_nm"]: float(r["E_time"]) for r in t},
                      "2600_endpoint_error_factor": float(next(r["endpoint_error_factor"] for r in h
                                                               if int(r["thickness_nm"]) == 2600))}
    comparisons = {}
    for complex_m, simple_m in (("M1", "M0"), ("M2", "M1")):
        a, b = metrics[complex_m], metrics[simple_m]
        comparisons[f"{complex_m}_vs_{simple_m}"] = {
            "time_reduction_fraction": 1 - a["E_time_macro"] / b["E_time_macro"],
            "time_films_improved": sum(a["per_film_time"][str(h)] < b["per_film_time"][str(h)]
                                        for h in hvalues),
            "thickness_change_fraction": a["E_h_macro"] / b["E_h_macro"] - 1,
            "passes_relative_thresholds": bool(a["E_time_macro"] <= 0.8 * b["E_time_macro"] and
                                               sum(a["per_film_time"][str(h)] < b["per_film_time"][str(h)]
                                                   for h in hvalues) >= 4 and
                                               a["E_h_macro"] <= 1.1 * b["E_h_macro"])}
    result = {"metrics": metrics, "comparisons": comparisons,
              "scope": "retrospective internal holdouts; not independent experiments",
              "design_thresholds": "relative 20% time, >=4 films, <=10% E_h increase; 2600 endpoint factor <=2 for quantitative outer extrapolation"}
    (OUT / "09_留出门槛核对.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), constrained_layout=True)
    x = np.arange(len(hvalues))
    colors = {"M0": "#596574", "M1": "#1d78b5", "M2": "#dc7c24"}
    for i, m in enumerate(models):
        axes[0].bar(x + (i - 1) * 0.24,
                    [metrics[m]["per_film_time"][str(h)] for h in hvalues],
                    width=0.23, color=colors[m], label=m)
        axes[1].bar(x + (i - 1) * 0.24,
                    [metrics[m]["per_film_h"][str(h)] for h in hvalues],
                    width=0.23, color=colors[m], label=m)
    for ax, title in zip(axes, ("Train 0–60 min, predict later times",
                                 "Leave one entire thickness out")):
        ax.set_title(title)
        ax.set_xticks(x, [str(h) for h in hvalues])
        ax.set_xlabel("Thickness (nm)")
        ax.set_ylabel("Retrospective asinh RMSE")
        ax.grid(axis="y", alpha=0.2)
        ax.legend(frameon=False)
    fig.savefig(OUT / "figures" / "时间与厚度留出逐膜误差.png", dpi=180)
    plt.close(fig)
    print(json.dumps(result["comparisons"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
