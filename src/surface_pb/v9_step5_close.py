"""Validate the predeclared MODEL_ONLY Step 5 stopping branch."""

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def rows(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    calibration = load("WP4/v9_第4步判定.json")
    model_map = load("WP5/v9_MODEL_ONLY_精简图判定.json")
    numerics = load("WP5/v9_MODEL_ONLY_十二状态数值复核.json")
    sensitivity = load("WP5/v9_MODEL_ONLY_八项敏感性判定.json")
    dielectric_check = load("WP5/v9_二维介电求解器核对.json")
    raw = rows("WP5/v9_MODEL_ONLY_八项敏感性.csv")
    if calibration["compatible_gamma_count"] != 0:
        raise AssertionError("wrong Step 5 branch")
    if model_map["reduced_map_fields"] != 34 or model_map["constant_charge_fields"] != 24:
        raise AssertionError("reduced map or fixed-charge control incomplete")
    if numerics["overall"] != "PASS" or numerics["states"] != 12:
        raise AssertionError("final twelve-state numerical check incomplete")
    if len(raw) != 144 or sensitivity["raw_rows"] != 144 or sensitivity["delta_rows"] != 126:
        raise AssertionError("frozen sensitivity matrix incomplete")
    if dielectric_check["status"] != "PASS":
        raise AssertionError("finite dielectric validation incomplete")
    if any(r["status"] != "VALID" for r in raw):
        raise AssertionError("unresolved sensitivity state")
    outcome = {
        "execution_status": "COMPLETE_PREDECLARED_MODEL_ONLY_BRANCH",
        "scientific_material_support": False,
        "calibration_state": calibration["calibration_state"],
        "calibration_E_min": calibration["calibration_E_min"],
        "full_A_grid_status": "NOT_RUN_PREDECLARED_GATE",
        "reduced_map_fields": model_map["reduced_map_fields"],
        "thick_N40_included_in_reduced_map": True,
        "fixed_charge_fields": model_map["constant_charge_fields"],
        "N4_N11_numerical_fields": model_map["N4_N11_extra_fields"],
        "N4_N11_abs_sigma_gain_C_m2": model_map["magnitude_difference_C_m2"],
        "N4_N11_five_times_error_C_m2": model_map["five_times_sum_error_C_m2"],
        "N4_N11_ten_percent_threshold_C_m2": model_map["ten_percent_sigma_N11_C_m2"],
        "N4_N11_model_only_enhancement_pass": model_map["model_only_enhancement_exceeds_both_thresholds"],
        "final_numerical_states": numerics["states"],
        "final_numerical_additional_PB_fields": numerics["new_PB_fields_total"]
                                                  + dielectric_check["additional_refinement_PB_fields"],
        "final_numerical_sigma_change_rel_max": numerics["max_epsilon_sigma_relative"],
        "final_numerical_zeta_change_mV_max": numerics["max_epsilon_zeta_mV"],
        "sensitivity_rows": sensitivity["raw_rows"],
        "sensitivity_equilibrium_PB_fields": 126,
        "sensitivity_H_only_observations": 18,
        "dielectric_max_sigma_change_relative": sensitivity["aggregates"]["PE_eps_2.3"]["max_abs_delta_sigma_relative"],
        "dielectric_max_zeta_change_mV": sensitivity["aggregates"]["PE_eps_2.3"]["max_abs_delta_zeta_mV"],
        "dielectric_crosses_frozen_risk_gate": sensitivity["dielectric_physics_unresolved_by_frozen_threshold"],
        "g_half_neutral_zeta_contrast_direction_changed": sensitivity["neutral_604_minus182_proxy_trends"]["g_x0.5"]["backed_zeta_direction_changed_from_baseline"],
        "g_half_neutral_zeta_contrast_below_1mV_numerical_tolerance":
            abs(sensitivity["neutral_604_minus182_proxy_trends"]["g_x0.5"]["backed_zeta_604_minus182_mV"]) < 1,
        "model_only_scope": "legacy example chemistry, proxy geometry, conditional NaCl; no material inference",
    }
    (ROOT / "WP5/v9_第5步验收.json").write_text(
        json.dumps(outcome, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"execution_status": outcome["execution_status"],
                      "scientific_material_support": outcome["scientific_material_support"],
                      "dielectric_crosses_frozen_risk_gate": outcome["dielectric_crosses_frozen_risk_gate"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
