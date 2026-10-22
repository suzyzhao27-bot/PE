"""Build WP4 comparison, salt budgets, convergence evidence, and figure."""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "WP4"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "effective_network_model"))
sys.path.insert(0, str(ROOT / "WP2" / "src"))
from transport import build_operator, solve  # noqa: E402
from v6_fit_wp4 import OUT, csv_write, config, load_films, predict  # noqa: E402


def main() -> None:
    films = load_films()
    models = ("M0", "M1", "M2")
    summaries = {m: json.loads((OUT / f"{m}_拟合摘要.json").read_text(encoding="utf-8")) for m in models}
    parameter_rows = []
    obs_rows = []
    comparison_rows = []
    budget_rows = []
    grid_rows = []
    predictions = {}

    for model in models:
        summary = summaries[model]
        p = summary["parameters"]
        parameter_rows.append({"model_id": model, "scenario_id": "main_r0_a25_rv1_dbref",
                               "best_start_id": summary["best_start_id"],
                               "loss_fit_eigen": summary["best_loss"],
                               "successful_starts": sum(x["success"] for x in summary["start_status"]),
                               "failed_starts": sum(not x["success"] for x in summary["start_status"]),
                               "boundary_hits": json.dumps(summary["boundary_hits"]),
                               **p, "raw_source_sha256_from_wp2": summary["raw_source_sha256_reused_from_wp2"]})
        for film in films:
            h = film.h_nm
            rows = list(csv.DictReader((OUT / f"{model}_观测预测.csv").open(encoding="utf-8-sig")))
            obs_rows += [r for r in rows if int(r["thickness_nm"]) == h]
            yobs = film.increment
            pred, info = predict(film, p, model, 50)
            predictions[(model, h)] = pred
            transformed = np.arcsinh(pred[1:]) - np.arcsinh(yobs[1:])
            early = film.times_min[1:] <= 60
            late = film.times_min[1:] > 60
            comparison_rows.append({"model_id": model, "thickness_nm": h,
                                    "loss_film": float(np.mean(transformed ** 2)),
                                    "rmse_asinh": float(np.sqrt(np.mean(transformed ** 2))),
                                    "rmse_early_asinh": float(np.sqrt(np.mean(transformed[early] ** 2))),
                                    "rmse_late_asinh": float(np.sqrt(np.mean(transformed[late] ** 2))),
                                    "mean_signed_residual_early": float(np.mean(transformed[early])),
                                    "mean_signed_residual_late": float(np.mean(transformed[late])),
                                    "endpoint_observed": yobs[-1], "endpoint_predicted": pred[-1],
                                    "endpoint_pred_over_obs": pred[-1] / yobs[-1],
                                    "mass_error_50": info["mass_error"],
                                    "minimum_concentration_50": info["minimum_concentration"],
                                    "numerical_status": "BDF success"})
            for n in (50, 100, 200):
                y, diagnostics = predict(film, p, model, n)
                numeric_loss = float(np.mean((np.arcsinh(y[1:]) - np.arcsinh(yobs[1:])) ** 2))
                grid_rows.append({"model_id": model, "thickness_nm": h, "n_cells": n,
                                  "loss_film": numeric_loss,
                                  "max_raw_difference_vs_50": float(np.max(np.abs(y - pred))),
                                  "max_scaled_difference_vs_50": float(np.max(np.abs(y - pred) / np.maximum(1, np.abs(pred)))),
                                  "mass_error": diagnostics["mass_error"],
                                  "minimum_concentration": diagnostics["minimum_concentration"]})
            if model != "M0":
                cfg = config(film, p)
                traj = solve(cfg, film.times_min * 60, n_cells=50, rtol=1e-7, atol=1e-10)
                op = traj.operator
                for j, t in enumerate(film.times_min):
                    state = traj.states[:, j]
                    n = op.n_cells
                    mobile = cfg.gamma_per_m * op.dx_m * np.sum(state[1:n + 1])
                    stored = cfg.gamma_per_m * op.dx_m * np.sum(state[n + 2:]) if op.has_storage else 0.0
                    budget_rows.append({"model_id": model, "thickness_nm": h, "time_min": t,
                                        "left_pool_mass_VR_C0": cfg.volume_ratio * state[0],
                                        "right_pool_mass_VR_C0": state[n + 1],
                                        "mobile_mass_VR_C0": mobile,
                                        "slow_mass_VR_C0": stored,
                                        "total_mass_VR_C0": op.mass(state),
                                        "gamma_h": cfg.gamma_per_m * cfg.h_m,
                                        "dm_over_db": p["dm_ratio"], "chi": p["chi"]})
            else:
                for j, t in enumerate(film.times_min):
                    right = pred[j] / p["q_conductivity"]
                    budget_rows.append({"model_id": model, "thickness_nm": h, "time_min": t,
                                        "left_pool_mass_VR_C0": 1 - right,
                                        "right_pool_mass_VR_C0": right,
                                        "mobile_mass_VR_C0": "not represented by M0",
                                        "slow_mass_VR_C0": "not represented by M0",
                                        "total_mass_VR_C0": 1,
                                        "gamma_h": p["gamma_per_m"] * h * 1e-9,
                                        "dm_over_db": p["dm_ratio"], "chi": 0})

    csv_write(OUT / "01_共享拟合参数表.csv", parameter_rows)
    csv_write(OUT / "02_全部观测与预测.csv", obs_rows)
    csv_write(OUT / "03_残差与模型比较.csv", comparison_rows)
    csv_write(OUT / "07_盐量时间账本.csv", budget_rows)
    csv_write(OUT / "08_拟合点数值加密.csv", grid_rows)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)
    palette = {"M0": "#5b6472", "M1": "#1c75bc", "M2": "#d77b24"}
    for ax, film in zip(axes.flat, films):
        ax.scatter(film.times_min, film.increment, c="black", s=30, label="Book1", zorder=4)
        for model in models:
            ax.plot(film.times_min, predictions[(model, film.h_nm)], "o-", ms=3,
                    lw=1.5, color=palette[model], label=model)
        ax.set_title(f"{film.h_nm} nm")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Conductivity increment (source unit)")
        ax.grid(alpha=0.2)
        if film.h_nm == 65:
            ax.legend(frameon=False)
    fig.suptitle("Frozen Book1 PE curves: shared-parameter training comparison", fontsize=13)
    fig.savefig(OUT / "figures" / "六厚度全曲线共同拟合.png", dpi=180)
    plt.close(fig)

    losses = {m: float(np.mean([r["loss_film"] for r in comparison_rows if r["model_id"] == m])) for m in models}
    report = {"loss_bdf_50": losses,
              "grid_max_scaled_50_to_200": {m: max(r["max_scaled_difference_vs_50"] for r in grid_rows
                                                if r["model_id"] == m and r["n_cells"] == 200) for m in models},
              "grid_max_mass_error": max(r["mass_error"] for r in grid_rows),
              "grid_minimum_state": min(r["minimum_concentration"] for r in grid_rows),
              "m0_m1_loss_reduction_fraction": 1 - losses["M1"] / losses["M0"],
              "m1_m2_loss_reduction_fraction": 1 - losses["M2"] / losses["M1"]}
    (OUT / "09_比较核对.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
