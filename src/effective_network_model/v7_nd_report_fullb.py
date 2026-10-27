"""Report completed 12-start ND full-B fit without claiming holdout success."""

from __future__ import annotations

import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_objective import H_VALUES, OUT, curves, residuals
from v7_nd_objective import BOUNDS, init_worker, predict_nd


def main():
    records = [json.loads(line) for line in
               (OUT / "v7_d1_fit_ND_fullB_runs.jsonl").read_text(
                   encoding="utf-8").splitlines()]
    if len(records) != 12 or sorted(r["start_index"] for r in records) != list(range(12)):
        raise RuntimeError("12 unique frozen ND full-B starts have not finished")
    good = [r for r in records if r["status"] == "solver_success" and
            r["best_L"] is not None]
    if not good:
        raise RuntimeError("No successful ND full-B start")
    stable = any(abs(a["best_L"] - b["best_L"]) /
                 min(a["best_L"], b["best_L"]) <= 0.01
                 for i, a in enumerate(good) for b in good[i + 1:])
    best = min(good, key=lambda r: r["best_L"])
    x = best["best_transformed"]
    alpha, q, beta, tau0, b = 10 ** x[0], 10 ** x[1], x[2], 10 ** x[3], x[4]
    data = curves()
    with ProcessPoolExecutor(max_workers=6, initializer=init_worker) as pool:
        predicted, diagnostics = predict_nd(pool, alpha, q, beta, tau0, b)
    residual = residuals(data, predicted)
    loss = float(residual @ residual)
    assert abs(loss - best["best_L"]) <= 1e-6 * max(1.0, best["best_L"])
    points, per_membrane = [], []
    for h in H_VALUES:
        baseline = data[h][0][1]
        obs, pred, rs = [], [], []
        for index, ((t, k, curve_id), y) in enumerate(
                zip(data[h], predicted[h], strict=True)):
            measured = k - baseline
            r = math.asinh(y) - math.asinh(measured)
            points.append({"model": "ND", "series": "B", "h_nm": h,
                           "curve_id": curve_id, "point_index": index, "time_min": t,
                           "conductivity_uS_cm": k, "baseline_uS_cm": baseline,
                           "observed_increment_uS_cm": measured,
                           "predicted_increment_uS_cm": y,
                           "prediction_seed_min_uS_cm": q * diagnostics[h][
                               "u_R_seed_min"][index],
                           "prediction_seed_max_uS_cm": q * diagnostics[h][
                               "u_R_seed_max"][index],
                           "asinh_residual": r, "in_training": index > 0})
            if index > 0:
                obs.append(measured); pred.append(y); rs.append(r)
        per_membrane.append({"h_nm": h, "n_train": len(obs),
                             "RMSE_original_uS_cm": float(np.sqrt(np.mean(
                                 np.square(np.array(pred) - np.array(obs))))),
                             "E_asinh": float(np.sqrt(np.mean(np.square(rs)))),
                             "endpoint_ratio_pred_over_obs": pred[-1] / obs[-1]
                             if obs[-1] != 0 else None,
                             "max_mass_relative_error": diagnostics[h][
                                 "max_mass_relative_error"],
                             "minimum_u": diagnostics[h]["minimum_u"]})
    with (OUT / "v7_ND_fullB_points.csv").open(
            "w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        writer.writeheader(); writer.writerows(points)
    endpoint_obs = {h: next(p for p in points if p["h_nm"] == h and
                            p["point_index"] == 6)["observed_increment_uS_cm"]
                    for h in (200, 350)}
    endpoint_pred = {h: next(p for p in points if p["h_nm"] == h and
                             p["point_index"] == 6)["predicted_increment_uS_cm"]
                     for h in (200, 350)}
    boundary = [name for name, value, (lo, hi) in zip(
        ("log10_alpha", "log10_Q", "beta", "log10_tau0_s", "b"), x, BOUNDS)
        if min(value - lo, hi - value) < 0.01 * (hi - lo)]
    summary = {"model": "ND",
               "status": "multi_start_stable" if stable else "local_optimum_unstable",
               "12_start_status_counts": {status: sum(r["status"] == status
                                                for r in records)
                                          for status in sorted(set(r["status"]
                                                                   for r in records))},
               "actual_objective_calls_total": sum(r["actual_objective_calls"]
                                                   for r in records),
               "best_start_index": best["start_index"],
               "best_transformed": x, "alpha_m_inv": alpha,
               "Q_uS_cm_per_unit_uR": q, "beta": beta, "tau0_s": tau0,
               "b": b, "boundary_within_one_percent": boundary,
               "L_recomputed": loss, "E_sqrt_L": math.sqrt(loss),
               "best_vs_worst_success_L_relative": (
                   max(r["best_L"] for r in good) - min(r["best_L"] for r in good)) /
                   min(r["best_L"] for r in good),
               "per_membrane": per_membrane,
               "worst_membrane_by_asinh_E": max(
                   per_membrane, key=lambda r: r["E_asinh"])["h_nm"],
               "endpoint_350_over_200_observed": endpoint_obs[350] / endpoint_obs[200],
               "endpoint_350_over_200_predicted": endpoint_pred[350] / endpoint_pred[200],
               "scope": "ND full B training only; T8 formal point and all holdouts pending"}
    (OUT / "v7_ND_fullB_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"ND full B: stable={stable}; E={summary['E_sqrt_L']:.6g}; "
          f"worst h={summary['worst_membrane_by_asinh_E']} nm")


if __name__ == "__main__":
    main()
