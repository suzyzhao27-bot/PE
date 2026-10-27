"""Compile the frozen WP4 comparison only after all N0/NS/ND reports exist."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP4"
MODELS = ("N0", "NS", "ND")
TASKS = ("fullB", "time", "film", "L")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    assert rows, path
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    needed = [OUT / f"v7_{model}_{task}_{kind}.{ext}"
              for model in MODELS for task in TASKS
              for kind, ext in (("summary", "json"), ("points", "csv"))]
    missing = [str(p.relative_to(ROOT)) for p in needed if not p.exists()]
    if missing:
        raise RuntimeError("Frozen reports incomplete: " + ", ".join(missing))
    summaries = {(m, t): json.loads((OUT / f"v7_{m}_{t}_summary.json").read_text(encoding="utf-8"))
                 for m in MODELS for t in TASKS}
    source = {(m, t): read_csv(OUT / f"v7_{m}_{t}_points.csv")
              for m in MODELS for t in TASKS}
    assert all(len(source[m, t]) == (70 if t == "L" else 42)
               for m in MODELS for t in TASKS)
    for task in TASKS:
        reference = None
        for model in MODELS:
            identity = [(p.get("hold_h_nm", p.get("h_nm")), p["curve_id"],
                         p["point_index"], p["time_min"], p["observed_increment_uS_cm"])
                        for p in source[model, task]]
            if reference is None:
                reference = identity
            else:
                assert identity == reference, (task, model, "observations differ")

    predictions = []
    for model in MODELS:
        for task in TASKS:
            for p in source[model, task]:
                index = int(p["point_index"])
                role = (p["role"] if "role" in p else
                        ("baseline" if index == 0 else
                         ("holdout" if task == "film" else "train")))
                if task == "fullB":
                    assert p["in_training"] == ("False" if index == 0 else "True")
                if task == "film":
                    assert p["in_training"] == "False"
                if task == "time":
                    assert role == ("baseline" if index == 0 else
                                    ("train" if float(p["time_min"]) <= 60 else "holdout"))
                if task == "L":
                    assert role == ("baseline" if index == 0 else
                                    ("early_train" if float(p["time_min"]) <= 240
                                     else "late_prediction"))
                observed = float(p["observed_increment_uS_cm"])
                predicted = float(p["predicted_increment_uS_cm"])
                asinh_r = math.asinh(predicted) - math.asinh(observed)
                assert math.isclose(asinh_r, float(p["asinh_residual"]),
                                    rel_tol=1e-10, abs_tol=1e-10)
                predictions.append({
                    "model": model, "task": task, "series": "L" if task == "L" else "B",
                    "hold_h_nm": p.get("hold_h_nm", "") if task == "film" else "",
                    "h_nm": p.get("h_nm", p.get("hold_h_nm")),
                    "curve_id": p["curve_id"], "point_index": index,
                    "time_min": p["time_min"], "role": role,
                    "fit_or_prediction": "baseline" if index == 0 else
                    ("fit" if role in ("train", "early_train") else "prediction"),
                    "conductivity_uS_cm": p.get("conductivity_uS_cm", ""),
                    "baseline_uS_cm": p.get("baseline_uS_cm", ""),
                    "observed_increment_uS_cm": observed,
                    "predicted_increment_uS_cm": predicted,
                    "prediction_seed_min_uS_cm": p["prediction_seed_min_uS_cm"],
                    "prediction_seed_max_uS_cm": p["prediction_seed_max_uS_cm"],
                    "original_residual_uS_cm": predicted-observed,
                    "asinh_residual": asinh_r,
                    "source_points": f"WP4/v7_{model}_{task}_points.csv",
                })
    assert len(predictions) == 588
    residuals = [p for p in predictions if p["point_index"] > 0]
    assert len(residuals) == 525

    params = []
    for model in MODELS:
        for task in TASKS:
            report = summaries[model, task]
            fits = report["folds"] if task == "film" else [report]
            for fit in fits:
                full = summaries[model, "fullB"]
                row = {
                    "model": model, "task": task,
                    "held_out_h_nm": fit.get("hold_h_nm", ""),
                    "status": fit["status"],
                    "best_start_index": fit["best_start_index"],
                    "alpha_m_inv": fit.get("alpha_m_inv", fit.get("alpha_L_m_inv", "")),
                    "Q_uS_cm_per_unit_uR": fit.get("Q_uS_cm_per_unit_uR",
                                                   fit.get("Q_L_uS_cm_per_unit_uR", "")),
                    "beta": fit.get("beta", full["beta"]),
                    "tau0_s": fit.get("tau0_s", full.get("tau0_s", "")) if model == "ND" else "",
                    "b": fit.get("b", full.get("b", "")) if model == "ND" else "",
                    "E_train_asinh": fit.get("E_sqrt_L", fit.get("E_train", fit.get("E_early", ""))),
                    "E_prediction_asinh": (fit.get("E_time", "") if task == "time" else
                                           fit.get("E_holdout", "") if task == "film" else
                                           fit.get("E_Llate", "") if task == "L" else ""),
                    "original_four_start_status_counts": json.dumps(
                        fit.get("initial_four_start_status_counts",
                                fit.get("four_start_status_counts",
                                        fit.get("12_start_status_counts", {}))), ensure_ascii=False),
                    "final_four_start_status_counts": json.dumps(
                        fit.get("four_start_status_counts",
                                fit.get("12_start_status_counts", {})), ensure_ascii=False),
                    "extension_status_counts": json.dumps(
                        fit.get("extension_status_counts", {}), ensure_ascii=False),
                    "boundary_within_one_percent": json.dumps(
                        fit.get("boundary_within_one_percent", []), ensure_ascii=False),
                    "source_summary": f"WP4/v7_{model}_{task}_summary.json",
                }
                params.append(row)
    assert len(params) == 27

    old = read_csv(OUT / "v6_metric_harmonization.csv")
    old_metrics = {(r["model"], r["split"]): float(r["v7_definition_six_film_RMS_E_asinh"])
                   for r in old}
    comparison = []
    for model in ("M0", "M1", "M2"):
        comparison.append({
            "model": model, "family": "historical_v6", "E_fullB_asinh": "",
            "E_time_asinh": old_metrics[model, "time"],
            "E_h_asinh": old_metrics[model, "film"], "E_Llate_asinh": "",
            "time_gain_vs_nested_pct": "", "film_gain_vs_nested_pct": "",
            "late_gain_vs_nested_pct": "", "time_ratio_to_v6_M1":
            old_metrics[model, "time"] / old_metrics["M1", "time"],
            "worst_holdout_h_nm": "", "2600_film_endpoint_pred_over_obs": "",
            "observed_200_over_350_at_120_min": "",
            "fullB_predicted_200_over_350_at_120_min": "",
            "fullB_original_RMSE_uS_cm": "", "film_original_RMSE_uS_cm": "",
            "numerical_status": "historical_technical_conditional" if model == "M2"
            else "historical_not_revalidated_on_v7_geometry",
            "D3_recommendation_gate": "historical_reference_only",
            "source": "WP4/v6_metric_harmonization.csv",
        })
    for model in MODELS:
        full, time, film, late = (summaries[model, t] for t in TASKS)
        baseline = "N0" if model == "NS" else "NS" if model == "ND" else None
        b_200 = next(p for p in source[model, "fullB"] if p["h_nm"] == "200"
                     and float(p["time_min"]) == 120)
        b_350 = next(p for p in source[model, "fullB"] if p["h_nm"] == "350"
                     and float(p["time_min"]) == 120)
        film_2600 = next(f for f in film["folds"] if f["hold_h_nm"] == 2600)
        b_res = [p["original_residual_uS_cm"] for p in predictions
                 if p["model"] == model and p["task"] == "fullB" and p["point_index"] > 0]
        h_res = [p["original_residual_uS_cm"] for p in predictions
                 if p["model"] == model and p["task"] == "film" and p["point_index"] > 0]
        valid = (all(summaries[model, t]["status"] == "multi_start_stable"
                     for t in ("fullB", "time", "L")) and
                 film["status"] == "all_folds_multi_start_stable")
        comparison.append({
            "model": model, "family": "v7_frozen_network",
            "E_fullB_asinh": full["E_sqrt_L"],
            "E_time_asinh": time["E_time"], "E_h_asinh": film["E_h"],
            "E_Llate_asinh": late["E_Llate"],
            "time_gain_vs_nested_pct": 100*(1-time["E_time"]/summaries[baseline, "time"]["E_time"])
            if baseline else "",
            "film_gain_vs_nested_pct": 100*(1-film["E_h"]/summaries[baseline, "film"]["E_h"])
            if baseline else "",
            "late_gain_vs_nested_pct": 100*(1-late["E_Llate"]/summaries[baseline, "L"]["E_Llate"])
            if baseline else "",
            "time_ratio_to_v6_M1": time["E_time"]/old_metrics["M1", "time"],
            "worst_holdout_h_nm": film["worst_holdout_membrane_by_asinh_E"],
            "2600_film_endpoint_pred_over_obs": film_2600["endpoint_ratio_pred_over_obs"],
            "observed_200_over_350_at_120_min":
            float(b_200["observed_increment_uS_cm"])/float(b_350["observed_increment_uS_cm"]),
            "fullB_predicted_200_over_350_at_120_min":
            float(b_200["predicted_increment_uS_cm"])/float(b_350["predicted_increment_uS_cm"]),
            "fullB_original_RMSE_uS_cm": math.sqrt(sum(r*r for r in b_res)/len(b_res)),
            "film_original_RMSE_uS_cm": math.sqrt(sum(r*r for r in h_res)/len(h_res)),
            "numerical_status": "all_reported_multistart_stable" if valid
            else "some_reported_fits_unstable",
            "D3_recommendation_gate": "nested_baseline" if model == "N0"
            else "fail_time_necessary_20_percent" if
            100*(1-time["E_time"]/summaries[baseline, "time"]["E_time"]) < 20
            else "requires_remaining_gates",
            "source": ";".join(f"WP4/v7_{model}_{t}_summary.json" for t in TASKS),
        })
    run_files = (
        "v7_d1_fit_fullB_runs.jsonl", "v7_d1_fit_time_runs.jsonl",
        "v7_d1_fit_L_runs.jsonl", "v7_d1_fit_film_runs.jsonl",
        "v7_d1_fit_ND_fullB_runs.jsonl", "v7_d1_fit_ND_time_runs.jsonl",
        "v7_d1_fit_ND_L_runs.jsonl", "v7_d1_fit_ND_film_runs.jsonl",
        "v7_d1_fit_film_extensions.jsonl",
    )
    run_rows = []
    seen_original = set()
    for filename in run_files:
        phase = "extension" if filename.endswith("extensions.jsonl") else "original"
        for lineno, line in enumerate((OUT / filename).read_text(encoding="utf-8").splitlines(), 1):
            run = json.loads(line)
            key = (run["model"], run["task"], run.get("hold_h_nm", ""), run["start_index"])
            if phase == "original":
                assert key not in seen_original, key
                seen_original.add(key)
            run_rows.append({
                "model": run["model"], "task": run["task"],
                "hold_h_nm": run.get("hold_h_nm", ""),
                "start_index": run["start_index"], "phase": phase,
                "status": run["status"], "initial_status": run.get("initial_status", ""),
                "best_L": run["best_L"],
                "actual_objective_calls": run.get("actual_objective_calls",
                                                  run.get("cumulative_actual_objective_calls", "")),
                "additional_objective_calls": run.get("additional_actual_objective_calls", ""),
                "message": run["message"],
                "source_jsonl": f"WP4/{filename}", "source_line": lineno,
            })
    assert len(seen_original) == 132 and len(run_rows) == 134
    write_csv(OUT / "parameters.csv", params)
    write_csv(OUT / "all_predictions.csv", predictions)
    write_csv(OUT / "residuals.csv", residuals)
    write_csv(OUT / "model_comparison.csv", comparison)
    write_csv(OUT / "optimization_statuses.csv", run_rows)
    print(f"WP4 compiled: parameters={len(params)} predictions={len(predictions)} "
          f"residuals={len(residuals)} comparison={len(comparison)} "
          f"optimizer_records={len(run_rows)}")


if __name__ == "__main__":
    main()
