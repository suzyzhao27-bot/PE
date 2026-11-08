"""Audit one frozen E3 four-start refit and its original-B comparison."""

from __future__ import annotations

import argparse
import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_objective import H_VALUES, curves, graphs
from v7_e2_inverse_start import OUT, init_worker, predict, residuals
from v7_e3_fit_start import data_for_task
from v7_network_transport import SPEC_PATH, fixed_bath_p, system


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("eta", type=float, choices=(0.01, 0.05))
    parser.add_argument("seed", type=int, choices=range(8101, 8107))
    args = parser.parse_args()
    key = (args.eta, args.seed)
    report_path = OUT / "v7_E3_task_reports.jsonl"
    if report_path.exists():
        for line in report_path.read_text(encoding="utf-8").splitlines():
            old = json.loads(line)
            if (old["eta"], old["seed"]) == key:
                raise RuntimeError("E3 task already reported")
    path = OUT / "v7_E3_fit_runs.jsonl"
    runs = [json.loads(line) for line in path.read_text(
        encoding="utf-8").splitlines()]
    group = [r for r in runs if (r["eta"], r["seed"]) == key]
    if len(group) != 4 or sorted(r["start_index"] for r in group) != list(range(4)):
        raise RuntimeError("Four unique frozen E3 starts have not finished")
    good = [r for r in group if r["status"] == "solver_success"
            and r["best_L"] is not None]
    if not good:
        raise RuntimeError("No successful E3 refit start")
    best = min(good, key=lambda r: r["best_L"])
    stable = any(r is not best and
                 abs(r["best_L"] - best["best_L"]) /
                 max(best["best_L"], np.finfo(float).tiny) <= 0.01
                 for r in good)
    if any(r["status"] in ("max_eval", "max_eval_or_solver_stop")
           and r["best_L"] is not None and r["best_L"] < best["best_L"]
           for r in group):
        raise RuntimeError("Lower-loss unfinished E3 start requires frozen extension review")
    with (OUT / "v7_E1_equivalent_selection.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        reference = next(r for r in csv.DictReader(stream)
                         if r["label"] == "fullB_best")
    x_ref = np.array([float(reference[k])
                      for k in ("log10_alpha", "log10_Q", "beta")])
    x_fit = np.array(best["best_transformed"], dtype=float)
    bounds_width = np.array([5.0, 6.0, 6.0])
    boundary = [name for i, name in
                enumerate(("log10_alpha", "log10_Q", "beta"))
                if min(x_fit[i], bounds_width[i] - x_fit[i]) <
                0.01 * bounds_width[i]]
    data = data_for_task(args.eta, args.seed)
    with ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                             initargs=(data,)) as pool:
        fitted, diagnostics = predict(
            pool, 10 ** x_fit[0], x_fit[2], 10 ** x_fit[1])
        original_fit, _ = predict(
            pool, 10 ** x_ref[0], x_ref[2], 10 ** x_ref[1])
    fit_residual = residuals(data, fitted)
    L = float(fit_residual @ fit_residual)
    assert math.isclose(L, best["best_L"], rel_tol=1e-6, abs_tol=1e-8)
    ref_residual = residuals(data, original_fit)
    L_ref = float(ref_residual @ ref_residual)
    max_mass = max(d[0] for d in diagnostics.values())
    minimum_u = min(d[1] for d in diagnostics.values())
    assert max_mass <= 1e-8 and minimum_u >= -1e-12
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    nets = graphs()
    original = curves()
    p_ratios = {}
    points = []
    original_error = []
    for h in H_VALUES:
        p_ref, p_fit = [], []
        for graph in nets[h]:
            cap_r, edge_r, area_r = system(
                graph, x_ref[2], 10 ** x_ref[0], spec)
            cap_f, edge_f, area_f = system(
                graph, x_fit[2], 10 ** x_fit[0], spec)
            p_ref.append(fixed_bath_p(cap_r, edge_r, area_r))
            p_fit.append(fixed_bath_p(cap_f, edge_f, area_f))
        p_ratios[str(h)] = float(np.mean(p_fit) / np.mean(p_ref))
        true_baseline = original[h][0][1]
        for i, ((time_min, raw, _), prediction) in enumerate(
                zip(original[h], fitted[h], strict=True)):
            original_increment = raw - true_baseline
            observed = float(data[h]["observed"][i])
            if i:
                original_error.append(math.asinh(prediction) -
                                      math.asinh(original_increment))
            points.append({
                "eta": args.eta, "seed": args.seed, "h_nm": h,
                "point_index": i, "time_min": time_min,
                "original_B_increment_uS_cm": original_increment,
                "perturbed_increment_uS_cm": observed,
                "fitted_increment_uS_cm": prediction,
            })
    report = {
        "model": "NS", "task": "E3_raw_reading_perturbation_four_start",
        "eta": args.eta, "seed": args.seed,
        "status": "multi_start_stable" if stable else "local_optimum_unstable",
        "start_status_counts": {status: sum(r["status"] == status for r in group)
                                for status in sorted({r["status"] for r in group})},
        "actual_objective_calls_total":
            sum(r["actual_objective_calls"] for r in group),
        "best_start_index": best["start_index"],
        "reference_transformed": x_ref.tolist(),
        "fit_transformed": x_fit.tolist(),
        "delta_transformed": (x_fit - x_ref).tolist(),
        "boundary_within_one_percent": boundary,
        "L_recomputed": L, "L_reference_on_perturbed": L_ref,
        "E_to_original_B": float(np.sqrt(np.mean(np.square(original_error)))),
        "P_fit_over_original_by_h": p_ratios,
        "max_mass_relative_error": max_mass, "minimum_u": minimum_u,
        "point_count": len(points),
    }
    point_path = OUT / "v7_E3_point_predictions.csv"
    with point_path.open("a", encoding="utf-8-sig" if not point_path.exists()
                         else "utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        if stream.tell() == 0:
            writer.writeheader()
        writer.writerows(points)
    with report_path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(report, ensure_ascii=False) + "\n")
    print(f"E3 {key}: {report['status']}; L={L:.8g}; "
          f"original_B_E={report['E_to_original_B']:.8g}; "
          f"points={len(points)}")


if __name__ == "__main__":
    main()
