"""Report ND independent time training and untouched t>60 min predictions."""

from __future__ import annotations

import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_objective import H_VALUES, OUT, curves
from v7_d1_time_fit_start import time_residuals
from v7_nd_objective import BOUNDS, init_worker, predict_nd


def main():
    records = [json.loads(line) for line in
               (OUT / "v7_d1_fit_ND_time_runs.jsonl").read_text(encoding="utf-8").splitlines()]
    selected = [r for r in records if r["model"] == "ND" and
                r["task"] == "time_B_train_le_60_min"]
    if len(selected) != 4 or sorted(r["start_index"] for r in selected) != list(range(4)):
        raise RuntimeError("Four unique frozen time starts have not finished")
    good = [r for r in selected if r["status"] == "solver_success" and
            r["best_L"] is not None]
    if not good:
        raise RuntimeError("No successful time start; retain the raw run records")
    stable = any(abs(a["best_L"] - b["best_L"]) /
                 min(a["best_L"], b["best_L"]) <= 0.01
                 for i, a in enumerate(good) for b in good[i + 1:])
    best = min(good, key=lambda r: r["best_L"])
    x = best["best_transformed"]
    alpha, q, beta, tau0, b = 10 ** x[0], 10 ** x[1], x[2], 10 ** x[3], x[4]
    data = curves()
    with ProcessPoolExecutor(max_workers=6, initializer=init_worker) as pool:
        predicted, diagnostics = predict_nd(pool, alpha, q, beta, tau0, b)
    train_r = time_residuals(data, predicted, train=True)
    hold_r = time_residuals(data, predicted, train=False)
    train_L = float(train_r @ train_r)
    assert abs(train_L - best["best_L"]) <= 1e-6 * max(1.0, best["best_L"])
    rows = []
    per_membrane = []
    for h in H_VALUES:
        points = data[h]
        baseline = points[0][1]
        train_obs, train_pred, train_res = [], [], []
        hold_obs, hold_pred, hold_res = [], [], []
        for index, ((t, k, curve_id), y) in enumerate(zip(points, predicted[h], strict=True)):
            measured = k - baseline
            r = math.asinh(y) - math.asinh(measured)
            role = "baseline" if index == 0 else ("train" if t <= 60 else "holdout")
            rows.append({"model": "ND", "series": "B", "h_nm": h,
                         "curve_id": curve_id, "point_index": index, "time_min": t,
                         "role": role, "conductivity_uS_cm": k,
                         "baseline_uS_cm": baseline, "observed_increment_uS_cm": measured,
                         "predicted_increment_uS_cm": y,
                         "prediction_seed_min_uS_cm": q * diagnostics[h]["u_R_seed_min"][index],
                         "prediction_seed_max_uS_cm": q * diagnostics[h]["u_R_seed_max"][index],
                         "asinh_residual": r})
            if role == "train":
                train_obs.append(measured); train_pred.append(y); train_res.append(r)
            elif role == "holdout":
                hold_obs.append(measured); hold_pred.append(y); hold_res.append(r)
        per_membrane.append({"h_nm": h, "n_train": len(train_obs),
                             "n_holdout": len(hold_obs),
                             "E_train_asinh": float(np.sqrt(np.mean(np.square(train_res)))),
                             "E_time_asinh": float(np.sqrt(np.mean(np.square(hold_res)))),
                             "RMSE_train_uS_cm": float(np.sqrt(np.mean(
                                 np.square(np.array(train_pred) - np.array(train_obs))))),
                             "RMSE_time_uS_cm": float(np.sqrt(np.mean(
                                 np.square(np.array(hold_pred) - np.array(hold_obs))))),
                             "endpoint_ratio_pred_over_obs": hold_pred[-1] / hold_obs[-1]
                             if hold_obs[-1] != 0 else None,
                             "max_mass_relative_error": diagnostics[h]["max_mass_relative_error"],
                             "minimum_u": diagnostics[h]["minimum_u"]})
    with (OUT / "v7_ND_time_points.csv").open(
            "w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    boundary = [name for name, value, (lo, hi) in zip(
        ("log10_alpha", "log10_Q", "beta", "log10_tau0_s", "b"), x, BOUNDS)
        if min(value - lo, hi - value) < 0.01 * (hi - lo)]
    summary = {"model": "ND", "status": "multi_start_stable" if stable else
               "local_optimum_unstable", "four_start_status_counts": {
                   status: sum(r["status"] == status for r in selected)
                   for status in sorted(set(r["status"] for r in selected))},
               "actual_objective_calls_total": sum(r["actual_objective_calls"]
                                                   for r in selected),
               "best_start_index": best["start_index"], "best_transformed": x,
               "alpha_m_inv": alpha, "Q_uS_cm_per_unit_uR": q, "beta": beta,
               "tau0_s": tau0, "b": b,
               "boundary_within_one_percent": boundary,
               "L_train_recomputed": train_L, "E_train": math.sqrt(train_L),
               "L_time": float(hold_r @ hold_r), "E_time": float(np.linalg.norm(hold_r)),
               "per_membrane": per_membrane,
               "worst_holdout_membrane_by_asinh_E": max(
                   per_membrane, key=lambda r: r["E_time_asinh"])["h_nm"],
               "scope": "B t<=60 min training, t>60 min time predictions; no full-data warm start"}
    (OUT / "v7_ND_time_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"ND time: stable={stable}; E_train={summary['E_train']:.6g}; "
          f"E_time={summary['E_time']:.6g}; "
          f"worst h={summary['worst_holdout_membrane_by_asinh_E']} nm")


if __name__ == "__main__":
    main()
