"""Eight frozen one-at-a-time sensitivities on nine actual proxy conditions."""

import csv
import json
import math
from pathlib import Path

from v9_geometry import family_geometry, make_geometry
from v9_observation import solve_velocity
from v9_pb_1d import Chemical
from v9_pb_2d import solve_equilibrium
from v9_pb_dielectric_2d import solve_dielectric

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP5/v9_MODEL_ONLY_八项敏感性.csv"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def save(rows):
    if not rows:
        return
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def geometry_change(spec, draw, kind, variation):
    baseline = family_geometry(spec, draw, kind)
    if variation in ("g_x0.5", "g_x2"):
        factor = 0.5 if variation == "g_x0.5" else 2.0
        gap = baseline.g_nm * factor
        N = 1 + math.floor((baseline.h_nm - baseline.d_nm) / (baseline.d_nm + gap) + 0.5)
        return make_geometry(kind, baseline.d_nm, gap, N, baseline.h_nm)
    if variation == "aligned_layers":
        return family_geometry(spec, draw, kind, aligned=True)
    return baseline


def solve_case(geometry, salt, pH, H, chemistry, variation):
    field = (solve_dielectric(geometry, salt, pH, chemistry, eps_PE_r=2.3)
             if variation == "PE_eps_2.3" else
             solve_equilibrium(geometry, salt, pH, chemical=chemistry))
    observation = solve_velocity(geometry, field, H) if H is not None else None
    return field, observation


def main():
    decision = json.loads((ROOT / "WP4/v9_第4步判定.json").read_text(encoding="utf-8"))
    if decision["compatible_gamma_count"] != 0:
        raise RuntimeError("this is the no-compatible MODEL_ONLY sensitivity branch")
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    old = spec["chemistry"]["legacy_fallback"]
    chemistry = Chemical(old["Gamma_nm2"], old["pI"], old["b"])
    observations = read_csv(ROOT / "WP1/v9_九组观测与条件盐度.csv")
    variations = ["baseline"] + spec["matrix"]["sensitivities"]
    rows = read_csv(OUT) if OUT.exists() else []
    for row in rows:
        N, h, d, period = int(row["N"]), float(row["h_nm"]), float(row["d_nm"]), float(row["period_nm"])
        row["vertical_spacing_nm"] = (h - d) / (N - 1) if N > 1 else "NA"
        row["model_water_area_fraction"] = 1 - N * d * d / (period * h)
    if rows:
        save(rows)
    existing = {(r["case_id"], r["kind"], r["variation"]) for r in rows}
    PB_fields = 0
    for source in observations:
        case_id = source["case_id"]
        draw = int(source["draw"])
        pH = float(source["pH_measured_mean"])
        base_salt = float(source["conditional_NaCl_mM"])
        base_H = float(source["H_full_gap_um_mean_assumed"])
        for kind in ("free", "backed"):
            base_geometry = family_geometry(spec, draw, kind)
            base_field = None
            for variation in variations:
                if variation.startswith("H_") and kind == "free":
                    continue
                key = (case_id, kind, variation)
                if key in existing:
                    continue
                geometry = geometry_change(spec, draw, kind, variation)
                salt = base_salt * (0.8 if variation == "cs_x0.8" else
                                    1.2 if variation == "cs_x1.2" else 1)
                H = None if kind == "free" else base_H * (0.9 if variation == "H_x0.9" else
                                                           1.1 if variation == "H_x1.1" else 1)
                row = {"case_id": case_id, "draw": draw, "pH_group": int(source["pH_group"]),
                       "kind": kind, "variation": variation,
                       "parameter_role": "MODEL_ONLY_LEGACY_REFERENCE_NOT_PE_FIT",
                       "salt_mM": salt, "pH_measured": pH,
                       "H_um": H if H is not None else "NA",
                       "N": geometry.N, "h_nm": geometry.h_nm, "d_nm": geometry.d_nm,
                       "g_nm": geometry.g_nm, "period_nm": geometry.period_nm,
                       "vertical_spacing_nm": geometry.vertical_spacing_nm if geometry.N > 1 else "NA",
                       "model_water_area_fraction": 1 - geometry.N * geometry.d_nm**2 /
                                                    (geometry.period_nm * geometry.h_nm),
                       "projected_open_fraction_proxy": geometry.g_nm / geometry.period_nm,
                       "wet_area_m2_per_m": geometry.wet_area_m2_per_m,
                       "aligned_layers": geometry.aligned}
                try:
                    if variation.startswith("H_"):
                        if base_field is None:
                            base_field, _ = solve_case(base_geometry, base_salt, pH, base_H,
                                                       chemistry, "baseline")
                            PB_fields += 1
                        field = base_field
                        observation = solve_velocity(base_geometry, field, H)
                    else:
                        field, observation = solve_case(geometry, salt, pH, H, chemistry, variation)
                        PB_fields += 1
                        if variation == "baseline":
                            base_field = field
                    row.update({"sigma_mean_C_m2": field.sigma_mean_C_m2,
                                "sigma_outer_C_m2": field.sigma_outer_C_m2,
                                "sigma_core_C_m2": field.sigma_core_C_m2 if field.sigma_core_C_m2 is not None else "NA",
                                "Q_proj_C_m2": field.Q_proj_C_m2,
                                "zeta_app_mV": observation.zeta_app_mV if observation else "NA",
                                "Gauss_absolute_C_m": field.Gauss_absolute_C_m,
                                "Gauss_relative": field.Gauss_relative if field.Gauss_relative is not None else "NA",
                                "max_species_M": field.max_species_M,
                                "min_species_M": field.min_species_M,
                                "physical_PB_warning": field.physical_PB_warning,
                                "layer_sigma_json": json.dumps(field.layer_sigma_C_m2, ensure_ascii=False),
                                "status": "VALID" if field.Gauss_relative is not None
                                                      and field.Gauss_relative <= 1e-6
                                                      and field.min_species_M >= 0
                                                      and not field.physical_PB_warning
                                          else "PHYSICAL_OR_GAUSS_RISK"})
                except Exception as exc:
                    row.update({"status": "NUMERICALLY_UNRESOLVED",
                                "error": f"{type(exc).__name__}: {exc}"})
                rows.append(row)
                existing.add(key)
                save(rows)
                print(f"{case_id} {kind} {variation}: {row['status']}", flush=True)
    expected = 9 * (2 * (1 + 6) + 2)
    if len(rows) != expected:
        raise AssertionError(f"expected {expected} sensitivity rows, got {len(rows)}")
    summary = {"rows": len(rows), "expected_rows": expected,
               "PB_fields_total": 126, "H_only_observations": 18,
               "PB_fields_computed_this_run": PB_fields,
               "valid_rows": sum(r["status"] == "VALID" for r in rows),
               "unresolved_rows": sum(r["status"] == "NUMERICALLY_UNRESOLVED" for r in rows),
               "max_Gauss_relative": max(float(r["Gauss_relative"]) for r in rows if r["status"] == "VALID"),
               "max_species_M": max(float(r["max_species_M"]) for r in rows if r["status"] == "VALID"),
               "step5_status": "SENSITIVITY_RAW_COMPLETE_PENDING_INTERPRETATION"}
    (ROOT / "WP5/v9_MODEL_ONLY_八项敏感性原始判定.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
