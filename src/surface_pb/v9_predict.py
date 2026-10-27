"""Lock 182 calibration witnesses, then evaluate frozen 357/604 holdouts once."""

import csv
import json
import math
from pathlib import Path

from v9_geometry import family_geometry
from v9_observation import solve_velocity
from v9_pb_1d import Chemical
from v9_pb_2d import solve_equilibrium

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP4"
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"


def read_jsonl(path):
    with path.open(encoding="utf-8") as stream:
        return [json.loads(line) for line in stream if line.strip()]


def write_csv(path, rows, columns):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    starts = read_jsonl(OUT / "v9_校准双起点.jsonl")
    gamma_grid = spec["chemistry"]["Gamma_nm2_grid"]
    assert len(starts) == 2 * len(gamma_grid)
    assert {(s["gamma_nm2"], s["start_index"]) for s in starts} == {
        (gamma, start) for gamma in gamma_grid for start in (1, 2)}
    assert all(s["actual_objective_calls"] <= 48 for s in starts)
    with (ROOT / "WP1/v9_九组观测与条件盐度.csv").open(encoding="utf-8-sig", newline="") as stream:
        obs = {r["case_id"]: r for r in csv.DictReader(stream)}
    by_gamma = []
    for gamma in gamma_grid:
        pair = [s for s in starts if s["gamma_nm2"] == gamma]
        valid = [s for s in pair if s["best_valid"]]
        best = min(valid, key=lambda s: s["best_valid"]["E_cal"]) if valid else None
        e_diff = abs(pair[0]["best_valid"]["E_cal"] - pair[1]["best_valid"]["E_cal"]) if len(valid) == 2 else None
        by_gamma.append({"gamma_nm2": gamma, "start1_status": pair[0]["status"],
                         "start2_status": pair[1]["status"],
                         "start1_E_cal": pair[0]["best_valid"]["E_cal"] if pair[0]["best_valid"] else "NA",
                         "start2_E_cal": pair[1]["best_valid"]["E_cal"] if pair[1]["best_valid"] else "NA",
                         "two_start_E_difference": e_diff if e_diff is not None else "NA",
                         "path_dependent_gt_0p1": bool(e_diff is not None and e_diff > 0.1),
                         "selected_start": best["start_index"] if best else "NA",
                         "pI": best["best_valid"]["pI"] if best else "NA",
                         "b": best["best_valid"]["b"] if best else "NA",
                         "E_cal": best["best_valid"]["E_cal"] if best else "NA",
                         "boundary_hit": best["boundary_hit"] if best else "NA",
                         "compatible": False})
    numeric = [r for r in by_gamma if r["E_cal"] != "NA"]
    e_min = min((r["E_cal"] for r in numeric), default=None)
    for row in by_gamma:
        row["compatible"] = bool(e_min is not None and row["E_cal"] != "NA"
                                 and row["E_cal"] <= 1 and row["E_cal"] <= e_min + 0.25)
    write_csv(OUT / "v9_校准七档汇总.csv", by_gamma, list(by_gamma[0]))
    compatible = [r for r in by_gamma if r["compatible"]]
    unresolved = any(s["status"] != "STABLE_LOCAL" for s in starts)
    ledger = {"calibration_E_min": e_min, "compatible_gamma_count": len(compatible),
              "calibration_state": ("COMPATIBLE_SET" if compatible else
                                    "OPTIMIZATION_UNRESOLVED_NO_COMPATIBLE_WITNESS" if unresolved else
                                    "PARAMETER_BOX_NO_COMPATIBLE_WITNESS"),
              "actual_objective_calls": sum(s["actual_objective_calls"] for s in starts),
              "maximum_PB_fields_from_objective_calls": 3 * sum(s["actual_objective_calls"] for s in starts),
              "prediction_set_chosen_without_holdouts": True}
    if not compatible:
        (OUT / "v9_第4步判定.json").write_text(json.dumps(ledger, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(ledger, ensure_ascii=False))
        return

    geom = {(draw, kind): family_geometry(spec, draw, kind)
            for draw in (182, 357, 604) for kind in ("backed", "free")}
    predictions = []
    for selected in compatible:
        chemical = Chemical(selected["gamma_nm2"], selected["pI"], selected["b"])
        for draw in (182, 357, 604):
            for pH_group in (3, 7, 10):
                row = obs[f"{draw}_pH{pH_group}"]
                for kind in ("backed", "free"):
                    if draw == 182 and kind == "backed":
                        continue  # the selected training call already contains these ζ values
                    entry = {"gamma_nm2": selected["gamma_nm2"], "pI": selected["pI"],
                             "b": selected["b"], "case_id": row["case_id"],
                             "draw": draw, "kind": kind, "role": row["role"],
                             "pH_measured": float(row["pH_measured_mean"]),
                             "conditional_NaCl_mM": float(row["conditional_NaCl_mM"]),
                             "zeta_obs_mV": float(row["zeta_reported_mean_mV"]) if kind == "backed" else "NA",
                             "comparison_scale_mV": float(row["comparison_scale_mV"]) if kind == "backed" else "NA"}
                    try:
                        field = solve_equilibrium(geom[(draw, kind)], entry["conditional_NaCl_mM"],
                                                  entry["pH_measured"], chemical=chemical)
                        entry.update({"sigma_mean_C_m2": field.sigma_mean_C_m2,
                                      "sigma_outer_C_m2": field.sigma_outer_C_m2,
                                      "sigma_core_C_m2": field.sigma_core_C_m2 if field.sigma_core_C_m2 is not None else "NA",
                                      "Q_proj_C_m2": field.Q_proj_C_m2,
                                      "Gauss_relative": field.Gauss_relative,
                                      "max_species_M": field.max_species_M,
                                      "physical_PB_warning": field.physical_PB_warning,
                                      "status": "VALID" if field.Gauss_relative <= 1e-6
                                                              and not field.physical_PB_warning else "PHYSICAL_OR_GAUSS_RISK"})
                        if kind == "backed":
                            observation = solve_velocity(geom[(draw, kind)], field,
                                                         float(row["H_full_gap_um_mean_assumed"]))
                            entry["zeta_pred_mV"] = observation.zeta_app_mV
                            entry["standardized_residual"] = ((observation.zeta_app_mV - entry["zeta_obs_mV"])
                                                               / entry["comparison_scale_mV"])
                        else:
                            entry["zeta_pred_mV"] = "NA"
                            entry["standardized_residual"] = "NA"
                    except Exception as exc:
                        entry["status"] = "NUMERICALLY_UNRESOLVED"
                        entry["error"] = f"{type(exc).__name__}: {exc}"
                    predictions.append(entry)
                    print(f"Γ={selected['gamma_nm2']:g} {row['case_id']} {kind} {entry['status']}", flush=True)
    columns = list(dict.fromkeys(key for row in predictions for key in row))
    write_csv(OUT / "v9_九组自由背衬预测.csv", predictions, columns)
    checks = []
    for selected in compatible:
        gamma = selected["gamma_nm2"]
        six = [r for r in predictions if r["gamma_nm2"] == gamma and r["kind"] == "backed"
               and r["draw"] != 182]
        valid = len(six) == 6 and all(r["status"] == "VALID" for r in six)
        E = math.sqrt(sum(r["standardized_residual"]**2 for r in six) / 6) if valid else None
        residual_gate = bool(valid and all(abs(r["standardized_residual"]) <= 2 for r in six))
        sign_gate = bool(valid and all((r["zeta_pred_mV"] > 0 if r["case_id"].endswith("pH3")
                                        else r["zeta_pred_mV"] < 0) for r in six))
        neutral604 = next((r for r in six if r["case_id"] == "604_pH7"), None)
        chosen_start = next(s for s in starts if s["gamma_nm2"] == gamma
                            and s["start_index"] == selected["selected_start"])
        train182_zeta = chosen_start["best_valid"]["zeta_pred_mV"][1]
        observed_difference = float(obs["604_pH7"]["zeta_reported_mean_mV"]) - float(obs["182_pH7"]["zeta_reported_mean_mV"])
        predicted_difference = neutral604["zeta_pred_mV"] - train182_zeta if valid else None
        difference_gate = bool(predicted_difference is not None
                               and abs(predicted_difference - observed_difference) <= 10)
        checks.append({"gamma_nm2": gamma, "E_cal": selected["E_cal"],
                       "E_pred": E if E is not None else "NA",
                       "max_abs_standardized_residual": max(abs(r["standardized_residual"]) for r in six) if valid else "NA",
                       "sign_gate": sign_gate, "residual_gate": residual_gate,
                       "observed_neutral_604_minus_182_mV": observed_difference,
                       "predicted_neutral_604_minus_182_mV": predicted_difference if predicted_difference is not None else "NA",
                       "neutral_difference_within_10mV": difference_gate,
                       "scientific_prediction_gate": bool(E is not None and E <= 1
                                                          and residual_gate and sign_gate and difference_gate)})
    write_csv(OUT / "v9_六组预测判定.csv", checks, list(checks[0]))
    ledger.update({"six_group_prediction_pass_gamma_count": sum(r["scientific_prediction_gate"] for r in checks),
                   "all_compatible_predictions_pass": all(r["scientific_prediction_gate"] for r in checks),
                   "predicted_fields": len(predictions)})
    (OUT / "v9_第4步判定.json").write_text(json.dumps(ledger, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(ledger, ensure_ascii=False))


if __name__ == "__main__":
    main()
