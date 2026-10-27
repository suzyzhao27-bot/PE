"""Report six frozen whole-membrane holdouts after all four starts per fold."""

from __future__ import annotations

import argparse
import csv
import json
import math

import numpy as np

from v7_d1_film_fit_start import film_residuals
from v7_d1_objective import H_VALUES, OUT, curves, graphs, predict
from v7_network_transport import SPEC_PATH


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("N0", "NS"))
    args = parser.parse_args()
    records = [json.loads(line) for line in
               (OUT / "v7_d1_fit_film_runs.jsonl").read_text(encoding="utf-8").splitlines()]
    extension_path = OUT / "v7_d1_fit_film_extensions.jsonl"
    extensions = ([json.loads(line) for line in
                   extension_path.read_text(encoding="utf-8").splitlines()]
                  if extension_path.exists() else [])
    data = curves()
    nets = graphs()
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    points, folds = [], []
    for hold_h in H_VALUES:
        selected = [r for r in records if r["model"] == args.model and
                    r["task"] == "B_whole_film_holdout" and r["hold_h_nm"] == hold_h]
        if len(selected) != 4 or sorted(r["start_index"] for r in selected) != list(range(4)):
            raise RuntimeError(f"Four unique starts not complete for hold_h={hold_h}")
        initial_status_counts = {status: sum(r["status"] == status for r in selected)
                                 for status in sorted({r["status"] for r in selected})}
        fold_extensions = [e for e in extensions if e["model"] == args.model and
                           e["hold_h_nm"] == hold_h]
        if len(fold_extensions) > 2 or len({e["start_index"] for e in fold_extensions}) != len(fold_extensions):
            raise RuntimeError(f"Invalid extension count for hold_h={hold_h}")
        unfinished = sorted((r for r in selected if
                             (r["status"] == "max_eval" or
                              (r["status"] == "max_eval_or_solver_stop" and
                               r.get("optimizer", {}).get("scipy_status") == 0)) and
                             r["best_L"] is not None),
                            key=lambda r: (r["best_L"], r["start_index"]))[:2]
        permitted = {r["start_index"] for r in unfinished}
        selected = [dict(r) for r in selected]
        for extension in fold_extensions:
            index = extension["start_index"]
            original = next(r for r in selected if r["start_index"] == index)
            if (index not in permitted or
                    extension["initial_actual_objective_calls"] != original["actual_objective_calls"] or
                    extension["initial_status"] != original["status"] or
                    extension["cumulative_actual_objective_calls"] > 400 or
                    extension["cumulative_actual_objective_calls"] !=
                    original["actual_objective_calls"] + extension["additional_actual_objective_calls"]):
                raise RuntimeError(f"Invalid extension provenance for hold_h={hold_h}, start={index}")
            original.update(status=extension["status"], best_L=extension["best_L"],
                            best_transformed=extension["best_transformed"],
                            actual_objective_calls=extension["cumulative_actual_objective_calls"])
        good = [r for r in selected if r["status"] == "solver_success" and
                r["best_L"] is not None]
        if not good:
            raise RuntimeError(f"No successful start for hold_h={hold_h}")
        stable = any(abs(a["best_L"] - b["best_L"]) /
                     min(a["best_L"], b["best_L"]) <= 0.01
                     for i, a in enumerate(good) for b in good[i + 1:])
        best = min(good, key=lambda r: r["best_L"])
        x = best["best_transformed"]
        alpha, q = 10 ** x[0], 10 ** x[1]
        beta = x[2] if args.model == "NS" else 0.0
        predicted, diagnostics = predict(data, nets, spec, alpha, beta, q)
        train_r = film_residuals(data, predicted, hold_h)
        train_L = float(train_r @ train_r)
        assert abs(train_L - best["best_L"]) <= 1e-6 * max(1.0, best["best_L"])
        baseline = data[hold_h][0][1]
        obs, pred, held_r = [], [], []
        for index, ((t, k, curve_id), y) in enumerate(
                zip(data[hold_h], predicted[hold_h], strict=True)):
            measured = k - baseline
            r = math.asinh(y) - math.asinh(measured)
            points.append({"model": args.model, "hold_h_nm": hold_h,
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
        bounds = [(0.0, 5.0), (0.0, 6.0)] + ([(0.0, 6.0)] if args.model == "NS" else [])
        boundary = [name for name, value, (lo, hi) in zip(
            ("log10_alpha", "log10_Q", "beta"), x, bounds)
            if min(value - lo, hi - value) < 0.01 * (hi - lo)]
        folds.append({"hold_h_nm": hold_h,
                      "status": "multi_start_stable" if stable else "local_optimum_unstable",
                      "initial_four_start_status_counts": initial_status_counts,
                      "extension_status_counts": {status: sum(e["status"] == status
                                                           for e in fold_extensions)
                                                  for status in sorted({e["status"]
                                                                        for e in fold_extensions})},
                      "four_start_status_counts": {status: sum(r["status"] == status
                                                       for r in selected)
                                                   for status in sorted(set(r["status"]
                                                                            for r in selected))},
                      "actual_objective_calls_total": sum(r["actual_objective_calls"]
                                                          for r in selected),
                      "best_start_index": best["start_index"],
                      "best_transformed": x, "alpha_m_inv": alpha,
                      "Q_uS_cm_per_unit_uR": q, "beta": beta,
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
    with (OUT / f"v7_{args.model}_film_points.csv").open(
            "w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        writer.writeheader(); writer.writerows(points)
    summary = {"model": args.model,
               "status": "all_folds_multi_start_stable" if all(
                   f["status"] == "multi_start_stable" for f in folds)
               else "some_folds_unstable",
               "E_h": float(np.sqrt(np.mean([f["E_holdout"] ** 2 for f in folds]))),
               "worst_holdout_membrane_by_asinh_E": max(
                   folds, key=lambda f: f["E_holdout"])["hold_h_nm"],
               "folds": folds,
               "scope": "six whole-membrane holdouts, each fitted on other five B curves"}
    (OUT / f"v7_{args.model}_film_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{args.model} film holdout: status={summary['status']}; "
          f"E_h={summary['E_h']:.6g}; "
          f"worst h={summary['worst_holdout_membrane_by_asinh_E']} nm")


if __name__ == "__main__":
    main()
