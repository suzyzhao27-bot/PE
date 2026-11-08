"""Report six frozen ND whole-membrane holdouts after four starts per fold."""

from __future__ import annotations

import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_film_fit_start import film_residuals
from v7_d1_objective import H_VALUES, OUT, curves
from v7_nd_objective import BOUNDS, init_worker, predict_nd


def main():
    records = [json.loads(line) for line in
               (OUT / "v7_d1_fit_ND_film_runs.jsonl").read_text(encoding="utf-8").splitlines()]
    data = curves()
    points, folds = [], []
    with ProcessPoolExecutor(max_workers=6, initializer=init_worker) as pool:
      for hold_h in H_VALUES:
        selected = [r for r in records if r["model"] == "ND" and
                    r["task"] == "B_whole_film_holdout" and r["hold_h_nm"] == hold_h]
        if len(selected) != 4 or sorted(r["start_index"] for r in selected) != list(range(4)):
            raise RuntimeError(f"Four unique starts not complete for hold_h={hold_h}")
        good = [r for r in selected if r["status"] == "solver_success" and
                r["best_L"] is not None]
        if not good:
            raise RuntimeError(f"No successful start for hold_h={hold_h}")
        best = min(good, key=lambda r: r["best_L"])
        any_pair_close = any(abs(a["best_L"] - b["best_L"]) /
                             min(a["best_L"], b["best_L"]) <= 0.01
                             for i, a in enumerate(good) for b in good[i + 1:])
        # The reported fit is the lowest-loss start; its own basin needs a
        # second successful start. A close pair in a worse basin cannot
        # certify the reported fit.
        stable = any(r is not best and
                     abs(r["best_L"] - best["best_L"]) / best["best_L"] <= 0.01
                     for r in good)
        x = best["best_transformed"]
        alpha, q, beta, tau0, b = 10 ** x[0], 10 ** x[1], x[2], 10 ** x[3], x[4]
        predicted, diagnostics = predict_nd(pool, alpha, q, beta, tau0, b)
        train_r = film_residuals(data, predicted, hold_h)
        train_L = float(train_r @ train_r)
        assert abs(train_L - best["best_L"]) <= 1e-6 * max(1.0, best["best_L"])
        baseline = data[hold_h][0][1]
        obs, pred, held_r = [], [], []
        for index, ((t, k, curve_id), y) in enumerate(
                zip(data[hold_h], predicted[hold_h], strict=True)):
            measured = k - baseline
            r = math.asinh(y) - math.asinh(measured)
            points.append({"model": "ND", "hold_h_nm": hold_h,
                           "curve_id": curve_id, "point_index": index, "time_min": t,
                           "observed_increment_uS_cm": measured,
                           "predicted_increment_uS_cm": y,
                           "prediction_seed_min_uS_cm": q * diagnostics[hold_h][
                               "u_R_seed_min"][index],
                           "prediction_seed_max_uS_cm": q * diagnostics[hold_h][
                               "u_R_seed_max"][index],
                           "asinh_residual": r, "in_training": False})
            if index > 0:
                obs.append(measured); pred.append(y); held_r.append(r)
        boundary = [name for name, value, (lo, hi) in zip(
            ("log10_alpha", "log10_Q", "beta", "log10_tau0_s", "b"), x, BOUNDS)
            if min(value - lo, hi - value) < 0.01 * (hi - lo)]
        folds.append({"hold_h_nm": hold_h,
                      "status": "multi_start_stable" if stable else "local_optimum_unstable",
                      "any_success_pair_within_one_percent": any_pair_close,
                      "selected_best_reproduced_within_one_percent": stable,
                      "four_start_status_counts": {status: sum(r["status"] == status
                                                       for r in selected)
                                                   for status in sorted(set(r["status"]
                                                                            for r in selected))},
                      "actual_objective_calls_total": sum(r["actual_objective_calls"]
                                                          for r in selected),
                      "best_start_index": best["start_index"],
                      "best_transformed": x, "alpha_m_inv": alpha,
                      "Q_uS_cm_per_unit_uR": q, "beta": beta,
                      "tau0_s": tau0, "b": b,
                      "boundary_within_one_percent": boundary,
                      "L_train_recomputed": train_L, "E_train": math.sqrt(train_L),
                      "E_holdout": float(np.sqrt(np.mean(np.square(held_r)))),
                      "RMSE_holdout_uS_cm": float(np.sqrt(np.mean(
                          np.square(np.array(pred) - np.array(obs))))),
                      "endpoint_ratio_pred_over_obs": pred[-1] / obs[-1]
                      if obs[-1] != 0 else None,
                      "max_mass_relative_error": diagnostics[hold_h][
                          "max_mass_relative_error"],
                      "minimum_u": diagnostics[hold_h]["minimum_u"]})
    with (OUT / "v7_ND_film_points.csv").open(
            "w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        writer.writeheader(); writer.writerows(points)
    summary = {"model": "ND",
               "status": "all_folds_multi_start_stable" if all(
                   f["status"] == "multi_start_stable" for f in folds)
               else "some_folds_unstable",
               "E_h": float(np.sqrt(np.mean([f["E_holdout"] ** 2 for f in folds]))),
               "worst_holdout_membrane_by_asinh_E": max(
                   folds, key=lambda f: f["E_holdout"])["hold_h_nm"],
               "folds": folds,
               "scope": "six whole-membrane holdouts, each fitted on other five B curves"}
    (OUT / "v7_ND_film_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"ND film holdout: status={summary['status']}; "
          f"E_h={summary['E_h']:.6g}; "
          f"worst h={summary['worst_holdout_membrane_by_asinh_E']} nm")


if __name__ == "__main__":
    main()
