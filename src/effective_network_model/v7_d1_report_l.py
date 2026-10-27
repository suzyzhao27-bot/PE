"""Report L early fit and late conditional-transfer predictions."""

from __future__ import annotations

import argparse
import csv
import json
import math

import numpy as np

from v7_d1_l_fit_start import H_L, curves_l, residuals_early
from v7_d1_objective import OUT, graphs
from v7_network_transport import SPEC_PATH, closed_baths, system


def predict_full(data, nets, spec, alpha, beta, q):
    predictions, diagnostics = {}, {}
    for h in H_L:
        times = np.array([t for t, _, _ in data[h]]) * 60.0
        runs = []
        for graph in nets[h]:
            capacities, edges, _ = system(graph, beta, alpha, spec)
            run = closed_baths(capacities, edges, times)
            runs.append(run)
        right = np.vstack([run["u_R"] for run in runs])
        predictions[h] = q * np.mean(right, axis=0)
        diagnostics[h] = {"max_mass_relative_error": max(
            run["mass_relative_error"] for run in runs),
            "minimum_u": min(run["minimum_u"] for run in runs),
            "seed_min": q * np.min(right, axis=0),
            "seed_max": q * np.max(right, axis=0)}
    return predictions, diagnostics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("N0", "NS"))
    args = parser.parse_args()
    records = [json.loads(line) for line in
               (OUT / "v7_d1_fit_L_runs.jsonl").read_text(encoding="utf-8").splitlines()]
    selected = [r for r in records if r["model"] == args.model and
                r["task"] == "L_early_le_240_min"]
    if len(selected) != 4 or sorted(r["start_index"] for r in selected) != list(range(4)):
        raise RuntimeError("Four unique frozen L starts have not finished")
    good = [r for r in selected if r["status"] == "solver_success" and
            r["best_L"] is not None]
    if not good:
        raise RuntimeError("No successful L start")
    stable = any(abs(a["best_L"] - b["best_L"]) /
                 min(a["best_L"], b["best_L"]) <= 0.01
                 for i, a in enumerate(good) for b in good[i + 1:])
    best = min(good, key=lambda r: r["best_L"])
    full_b = json.loads((OUT / f"v7_{args.model}_fullB_summary.json").read_text(
        encoding="utf-8"))
    assert full_b["status"] == "multi_start_stable"
    beta = float(full_b["beta"])
    assert all(abs(r["fixed_beta_from_full_B"] - beta) < 1e-12 for r in selected)
    x = best["best_transformed"]
    alpha, q = 10 ** x[0], 10 ** x[1]
    data = curves_l()
    nets = graphs()
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    predictions, diagnostics = predict_full(data, nets, spec, alpha, beta, q)
    early_predictions = {h: predictions[h][[i for i, (t, _, _) in
                                               enumerate(data[h]) if t <= 240]]
                         for h in H_L}
    train_r = residuals_early(data, early_predictions)
    train_L = float(train_r @ train_r)
    assert abs(train_L - best["best_L"]) <= 1e-6 * max(1.0, best["best_L"])
    points, per_membrane = [], []
    for h in H_L:
        baseline = data[h][0][1]
        final_obs = data[h][-1][1] - baseline
        final_pred = predictions[h][-1]
        early_r, late_r, late_obs, late_pred = [], [], [], []
        for index, ((t, k, curve_id), y) in enumerate(
                zip(data[h], predictions[h], strict=True)):
            observed = k - baseline
            r = math.asinh(y) - math.asinh(observed)
            role = "baseline" if index == 0 else ("early_train" if t <= 240
                                                  else "late_prediction")
            points.append({"model": args.model, "h_nm": h, "curve_id": curve_id,
                           "point_index": index, "time_min": t, "role": role,
                           "conductivity_uS_cm": k, "baseline_uS_cm": baseline,
                           "observed_increment_uS_cm": observed,
                           "predicted_increment_uS_cm": y,
                           "prediction_seed_min_uS_cm": diagnostics[h]["seed_min"][index],
                           "prediction_seed_max_uS_cm": diagnostics[h]["seed_max"][index],
                           "asinh_residual": r,
                           "observed_fraction_of_observed_endpoint_descriptive":
                           observed / final_obs if final_obs != 0 else None,
                           "predicted_fraction_of_predicted_endpoint_descriptive":
                           y / final_pred if final_pred != 0 else None})
            if role == "early_train":
                early_r.append(r)
            elif role == "late_prediction":
                late_r.append(r); late_obs.append(observed); late_pred.append(y)
        per_membrane.append({"h_nm": h, "n_early_train": len(early_r),
                             "n_late_prediction": len(late_r),
                             "E_early": float(np.sqrt(np.mean(np.square(early_r)))),
                             "E_late": float(np.sqrt(np.mean(np.square(late_r)))),
                             "RMSE_late_uS_cm": float(np.sqrt(np.mean(np.square(
                                 np.array(late_pred) - np.array(late_obs))))),
                             "endpoint_ratio_pred_over_obs": final_pred / final_obs
                             if final_obs != 0 else None,
                             "max_mass_relative_error": diagnostics[h][
                                 "max_mass_relative_error"],
                             "minimum_u": diagnostics[h]["minimum_u"]})
    with (OUT / f"v7_{args.model}_L_points.csv").open(
            "w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        writer.writeheader(); writer.writerows(points)
    boundary = [name for name, value, (lo, hi) in zip(
        ("log10_alpha_L", "log10_Q_L"), x, ((0.0, 5.0), (0.0, 6.0)))
        if min(value - lo, hi - value) < 0.01 * (hi - lo)]
    summary = {"model": args.model,
               "status": "multi_start_stable" if stable else "local_optimum_unstable",
               "four_start_status_counts": {status: sum(r["status"] == status
                                                for r in selected)
                                            for status in sorted(set(r["status"]
                                                                     for r in selected))},
               "actual_objective_calls_total": sum(r["actual_objective_calls"]
                                                   for r in selected),
               "best_start_index": best["start_index"],
               "best_transformed": x, "alpha_L_m_inv": alpha,
               "Q_L_uS_cm_per_unit_uR": q, "fixed_beta_from_full_B": beta,
               "boundary_within_one_percent": boundary,
               "L_early_recomputed": train_L, "E_early": math.sqrt(train_L),
               "E_Llate": float(np.sqrt(np.mean([r["E_late"] ** 2
                                                for r in per_membrane]))),
               "per_membrane": per_membrane,
               "scope": "L early<=240 min fitted; late>240 min conditional transfer; endpoint-normalized shapes descriptive only"}
    (OUT / f"v7_{args.model}_L_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{args.model} L: stable={stable}; E_early={summary['E_early']:.6g}; "
          f"E_Llate={summary['E_Llate']:.6g}")


if __name__ == "__main__":
    main()
