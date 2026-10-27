"""Validate one NS E1 profile point after its three frozen starts complete."""

from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor

from v7_d1_objective import ROOT, curves, init_forward_worker, predict_parallel, residuals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("parameter_index", type=int, choices=range(3))
    parser.add_argument("position_index", type=int, choices=range(9))
    args = parser.parse_args()
    path = ROOT / "WP5/v7_E1_profile_runs.jsonl"
    raw = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    selected = [r for r in raw if (r["parameter_index"], r["position_index"]) ==
                (args.parameter_index, args.position_index)]
    if len(selected) != 3 or sorted(r["start_index"] for r in selected) != [0, 1, 2]:
        raise RuntimeError("All three unique frozen profile starts must complete first")
    report_path = ROOT / "WP5/v7_E1_profile_point_reports.jsonl"
    if report_path.exists():
        for line in report_path.read_text(encoding="utf-8").splitlines():
            old = json.loads(line)
            if (old["parameter_index"], old["position_index"]) == (
                    args.parameter_index, args.position_index):
                raise RuntimeError("Profile point already reported")
    lower, upper = (0., 0., 0.), (5., 6., 6.)
    fixed = lower[args.parameter_index] + (upper[args.parameter_index] - lower[args.parameter_index]) * args.position_index / 8
    assert all(r["model"] == "NS" and
               math.isclose(r["fixed_value_transformed"], fixed, abs_tol=1e-12) and
               r["actual_objective_calls"] <= 200 for r in selected)
    original_counts = {status: sum(r["status"] == status for r in selected)
                       for status in sorted({r["status"] for r in selected})}
    original_unfinished = sorted((r for r in selected if
                                  (r["status"] == "max_eval" or
                                   (r["status"] == "max_eval_or_solver_stop" and
                                    r.get("optimizer", {}).get("scipy_status") == 0)) and
                                  r["best_L"] is not None and math.isfinite(r["best_L"])),
                                 key=lambda r: (r["best_L"], r["start_index"]))[:2]
    eligible = {r["start_index"] for r in original_unfinished}
    extension_path = ROOT / "WP5/v7_E1_profile_extensions.jsonl"
    extensions = ([json.loads(line) for line in extension_path.read_text(encoding="utf-8").splitlines()]
                  if extension_path.exists() else [])
    extensions = [e for e in extensions if (e["parameter_index"], e["position_index"]) ==
                  (args.parameter_index, args.position_index)]
    if len(extensions) > 2 or len({e["start_index"] for e in extensions}) != len(extensions):
        raise RuntimeError("Invalid profile extension count")
    selected = [dict(r) for r in selected]
    for extension in extensions:
        index = extension["start_index"]
        if index not in eligible:
            raise RuntimeError("Extension start is not among two lowest original unfinished starts")
        original = next(r for r in selected if r["start_index"] == index)
        if (extension["initial_actual_objective_calls"] != original["actual_objective_calls"] or
                extension["initial_status"] != original["status"] or
                extension["cumulative_actual_objective_calls"] > 400 or
                extension["cumulative_actual_objective_calls"] !=
                original["actual_objective_calls"] + extension["additional_actual_objective_calls"] or
                extension["best_L"] > original["best_L"] + 1e-8):
            raise RuntimeError("Invalid profile extension provenance")
        original.update(status=extension["status"], best_L=extension["best_L"],
                        best_transformed=extension["best_transformed"],
                        actual_objective_calls=extension["cumulative_actual_objective_calls"])
    good = [r for r in selected if r["status"] == "solver_success" and
            r["best_L"] is not None and math.isfinite(r["best_L"]) and
            r["best_transformed"] is not None]
    best = min(good, key=lambda r: r["best_L"]) if good else None
    stable = bool(best and any(r is not best and
                               abs(r["best_L"]-best["best_L"])/best["best_L"] <= 0.01
                               for r in good))
    unfinished = [r for r in selected if r["best_L"] is not None and
                  (r["status"] == "max_eval" or
                   (r["status"] == "max_eval_or_solver_stop" and
                    r.get("optimizer", {}).get("scipy_status") == 0))]
    lower_unfinished = bool(best and unfinished and
                            min(r["best_L"] for r in unfinished) < best["best_L"])
    pending = eligible - {e["start_index"] for e in extensions}
    if pending and (not stable or lower_unfinished):
        raise RuntimeError("Check the two lowest-loss unfinished starts for the sole 400-call extension before reporting")
    stable = stable and not lower_unfinished
    recomputed = None
    mass_error = None
    minimum_u = None
    if best:
        x = best["best_transformed"]
        assert len(x) == 3 and math.isclose(x[args.parameter_index], fixed, abs_tol=1e-12)
        with ProcessPoolExecutor(max_workers=6, initializer=init_forward_worker) as pool:
            predicted, diagnostics = predict_parallel(pool, 10**x[0], x[2], 10**x[1])
        r = residuals(curves(), predicted)
        recomputed = float(r @ r)
        mass_error = max(d["max_mass_relative_error"] for d in diagnostics.values())
        minimum_u = min(d["minimum_u"] for d in diagnostics.values())
        assert math.isclose(recomputed, best["best_L"], rel_tol=1e-6, abs_tol=1e-8)
        assert mass_error <= 1e-8 and minimum_u >= -1e-12
    report = {
        "model": "NS", "task": "E1_fixed_parameter_profile_point",
        "parameter_index": args.parameter_index,
        "parameter": selected[0]["parameter"],
        "position_index": args.position_index,
        "fixed_value_transformed": fixed,
        "status": "multi_start_stable" if stable else
                  "local_optimum_unstable" if best else "no_successful_start",
        "start_status_counts": {status: sum(r["status"] == status for r in selected)
                                for status in sorted({r["status"] for r in selected})},
        "original_start_status_counts": original_counts,
        "extension_status_counts": {status: sum(e["status"] == status for e in extensions)
                                    for status in sorted({e["status"] for e in extensions})},
        "best_start_index": best["start_index"] if best else None,
        "best_transformed": best["best_transformed"] if best else None,
        "best_L": best["best_L"] if best else None,
        "L_recomputed": recomputed,
        "max_mass_relative_error": mass_error,
        "minimum_u": minimum_u,
        "actual_objective_calls_total": sum(r["actual_objective_calls"] for r in selected),
    }
    with report_path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(report, ensure_ascii=False) + "\n")
    print(f"NS E1 {report['parameter']} position {args.position_index}: "
          f"{report['status']}; L={report['best_L']}; calls={report['actual_objective_calls_total']}")


if __name__ == "__main__":
    main()
