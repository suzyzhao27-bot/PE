"""Audit one frozen E2 four-start inversion, curves, and fixed-bath P."""

from __future__ import annotations

import argparse
import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_objective import H_VALUES, ROOT, graphs
from v7_e2_inverse_start import OUT, data_for_task, init_worker, predict, residuals
from v7_network_transport import SPEC_PATH, fixed_bath_p, system


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("truth_label", choices=("fullB_best", "equivalent_1"))
    parser.add_argument("noise_eta", type=float, choices=(0.0, 0.01))
    parser.add_argument("sampling", choices=("sparse_actual", "dense_1min"))
    parser.add_argument("startset_index", type=int, choices=(0, 1))
    args = parser.parse_args()
    key = (args.truth_label, args.noise_eta, args.sampling,
           args.startset_index)
    report_path = OUT / "v7_E2_task_reports.jsonl"
    if report_path.exists():
        for line in report_path.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            if (r["truth_label"], r["noise_eta"], r["sampling"],
                    r["startset_index"]) == key:
                raise RuntimeError("E2 task already reported")
    runs = [json.loads(line) for line in
            (OUT / "v7_E2_inverse_runs.jsonl").read_text(
                encoding="utf-8").splitlines()]
    group = [r for r in runs if
             (r["truth_label"], r["noise_eta"], r["sampling"],
              r["startset_index"]) == key]
    if len(group) != 4 or sorted(r["start_index"] for r in group) != list(range(4)):
        raise RuntimeError("Four unique frozen E2 starts have not finished")
    good = [r for r in group if r["status"] == "solver_success"
            and r["best_L"] is not None]
    if not good:
        raise RuntimeError("No successful E2 inverse start")
    best = min(good, key=lambda r: r["best_L"])
    stable = any(r is not best and
                 abs(r["best_L"] - best["best_L"]) /
                 max(best["best_L"], np.finfo(float).tiny) <= 0.01
                 for r in good)
    lower_unfinished = any(
        r["status"] in ("max_eval", "max_eval_or_solver_stop")
        and r["best_L"] is not None and r["best_L"] < best["best_L"]
        for r in group)
    if lower_unfinished:
        raise RuntimeError("Lower-loss unfinished E2 start requires frozen extension review")
    with (OUT / "v7_E1_equivalent_selection.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        truth = next(r for r in csv.DictReader(stream)
                     if r["label"] == args.truth_label)
    x_truth = np.array([float(truth[k])
                        for k in ("log10_alpha", "log10_Q", "beta")])
    x_fit = np.array(best["best_transformed"], dtype=float)
    data = data_for_task(args.truth_label, args.noise_eta, args.sampling)
    with ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                             initargs=(data,)) as pool:
        predicted, diagnostics = predict(
            pool, 10 ** x_fit[0], x_fit[2], 10 ** x_fit[1])
    r = residuals(data, predicted)
    L = float(r @ r)
    assert math.isclose(L, best["best_L"], rel_tol=1e-6, abs_tol=1e-8)
    max_mass = max(d[0] for d in diagnostics.values())
    minimum_u = min(d[1] for d in diagnostics.values())
    assert max_mass <= 1e-8 and minimum_u >= -1e-12
    with (OUT / "v7_E2_synthetic_observations.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        synthetic = [row for row in csv.DictReader(stream)
                     if (row["truth_label"], float(row["noise_eta"]),
                         row["sampling"]) ==
                     (args.truth_label, args.noise_eta, args.sampling)]
    points, curve_true_errors, p_ratios = [], [], {}
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    nets = graphs()
    for h in H_VALUES:
        rows = sorted((row for row in synthetic if int(row["h_nm"]) == h),
                      key=lambda row: float(row["time_min"]))
        assert len(rows) == len(predicted[h])
        true_baseline = float(rows[0]["noiseless_raw_kappa_uS_cm"])
        for i, (row, y) in enumerate(zip(rows, predicted[h], strict=True)):
            true_increment = (float(row["noiseless_raw_kappa_uS_cm"]) -
                              true_baseline)
            observed = float(row["increment_uS_cm"])
            curve_true_errors.append(math.asinh(y) - math.asinh(true_increment))
            points.append({
                "truth_label": args.truth_label, "noise_eta": args.noise_eta,
                "sampling": args.sampling,
                "startset_index": args.startset_index,
                "h_nm": h, "point_index": i,
                "time_min": float(row["time_min"]),
                "observed_increment_uS_cm": observed,
                "noiseless_true_increment_uS_cm": true_increment,
                "fitted_increment_uS_cm": y,
            })
        p_true, p_fit = [], []
        for graph in nets[h]:
            cap_t, edge_t, area_t = system(
                graph, x_truth[2], 10 ** x_truth[0], spec)
            cap_f, edge_f, area_f = system(
                graph, x_fit[2], 10 ** x_fit[0], spec)
            p_true.append(fixed_bath_p(cap_t, edge_t, area_t))
            p_fit.append(fixed_bath_p(cap_f, edge_f, area_f))
        mean_true = float(np.mean(p_true))
        mean_fit = float(np.mean(p_fit))
        p_ratios[str(h)] = (mean_fit / mean_true if mean_true > 0 else None)
    point_path = OUT / "v7_E2_inverse_point_predictions.csv"
    with point_path.open("a", encoding="utf-8-sig" if not point_path.exists()
                         else "utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        if stream.tell() == 0:
            writer.writeheader()
        writer.writerows(points)
    report = {
        "model": "NS", "task": "E2_synthetic_inverse_four_start",
        "truth_label": args.truth_label, "noise_eta": args.noise_eta,
        "sampling": args.sampling, "startset_index": args.startset_index,
        "status": "multi_start_stable" if stable else "local_optimum_unstable",
        "start_status_counts": {status: sum(r["status"] == status for r in group)
                                for status in sorted({r["status"] for r in group})},
        "actual_objective_calls_total": sum(r["actual_objective_calls"]
                                            for r in group),
        "best_start_index": best["start_index"],
        "truth_transformed": x_truth.tolist(),
        "fit_transformed": x_fit.tolist(),
        "absolute_transformed_errors": np.abs(x_fit - x_truth).tolist(),
        "distance_transformed": float(np.linalg.norm(x_fit - x_truth)),
        "L_recomputed": L,
        "E_to_noisy_observation": math.sqrt(L),
        "E_to_noiseless_truth_all_points":
            float(np.sqrt(np.mean(np.square(curve_true_errors)))),
        "P_fit_over_truth_by_h": p_ratios,
        "max_mass_relative_error": max_mass,
        "minimum_u": minimum_u,
        "point_count": len(points),
    }
    with report_path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(report, ensure_ascii=False) + "\n")
    print(f"E2 {key}: {report['status']}; L={L:.8g}; "
          f"distance={report['distance_transformed']:.8g}; "
          f"points={len(points)}")


if __name__ == "__main__":
    main()
