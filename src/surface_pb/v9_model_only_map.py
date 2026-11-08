"""Predeclared no-compatible-parameter reduced map and constant-charge control."""

import csv
import json
import math
from pathlib import Path

from v9_geometry import pure_geometry
from v9_observation import solve_velocity
from v9_pb_1d import Chemical
from v9_pb_2d import solve_equilibrium

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP5"
MAP = OUT / "v9_MODEL_ONLY_精简厚度图.csv"
CONTROL = OUT / "v9_MODEL_ONLY_固定电荷对照.csv"
ERROR = OUT / "v9_MODEL_ONLY_N4_N11数值误差.csv"


def read_csv(path):
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def append_csv(path, row):
    with path.open("a", encoding="utf-8-sig" if not path.exists() else "utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        if stream.tell() == 0:
            writer.writeheader()
        writer.writerow(row)


def calculate(spec, chemistry, N, salt, kind, fixed_sigma=None, mesh_step_nm=None,
              outer_factor=1.0):
    geometry = pure_geometry(spec, N, kind)
    field = solve_equilibrium(geometry, salt, 7.0,
                              chemical=chemistry if fixed_sigma is None else None,
                              fixed_sigma=fixed_sigma,
                              mesh_step_nm=mesh_step_nm, outer_factor=outer_factor)
    obs = solve_velocity(geometry, field, 100.0) if kind == "backed" else None
    return {"N": N, "h_nm": geometry.h_nm, "d_nm": geometry.d_nm,
            "g_nm": geometry.g_nm, "salt_mM": salt, "pH": 7.0, "kind": kind,
            "sigma_mean_C_m2": field.sigma_mean_C_m2,
            "sigma_outer_C_m2": field.sigma_outer_C_m2,
            "sigma_core_C_m2": field.sigma_core_C_m2 if field.sigma_core_C_m2 is not None else "NA",
            "Q_proj_C_m2": field.Q_proj_C_m2,
            "zeta_app_mV": obs.zeta_app_mV if obs else "NA",
            "Gauss_absolute_C_m": field.Gauss_absolute_C_m,
            "Gauss_relative": field.Gauss_relative if field.Gauss_relative is not None else "NA",
            "max_species_M": field.max_species_M,
            "physical_PB_warning": field.physical_PB_warning,
            "layer_sigma_json": json.dumps(field.layer_sigma_C_m2, ensure_ascii=False),
            "mesh_step_nm": field.mesh.target_nm, "L_nm": field.mesh.L_nm,
            "status": "VALID" if field.Gauss_relative is not None
                                 and field.Gauss_relative <= 1e-6 and not field.physical_PB_warning
                      else "PHYSICAL_OR_GAUSS_RISK"}


def main():
    decision = json.loads((ROOT / "WP4/v9_第4步判定.json").read_text(encoding="utf-8"))
    if decision["compatible_gamma_count"] != 0:
        raise RuntimeError("MODEL_ONLY fallback applies only with no compatible chemistry")
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    old = spec["chemistry"]["legacy_fallback"]
    chemistry = Chemical(old["Gamma_nm2"], old["pI"], old["b"])
    matrix = spec["matrix"]
    salts = matrix["pure_h_salt_mM"]
    map_cases = [(N, salt, kind) for salt in salts for N in (1, 4, 11, 20)
                 for kind in ("free", "backed")]
    map_cases += [(40, 1, kind) for kind in ("free", "backed")]
    existing = {(int(r["N"]), float(r["salt_mM"]), r["kind"]) for r in read_csv(MAP)}
    for N, salt, kind in map_cases:
        if (N, float(salt), kind) in existing:
            continue
        row = calculate(spec, chemistry, N, salt, kind)
        row["parameter_role"] = "MODEL_ONLY_LEGACY_REFERENCE_NOT_PE_FIT"
        append_csv(MAP, row)
        print(f"map N{N} c={salt:g} {kind} {row['status']}", flush=True)
    mapped = read_csv(MAP)
    assert len(mapped) == 34

    # N4/N11: estimate the declared difference error from half-step and doubled L.
    existing_error = {(int(r["N"]), r["variation"]) for r in read_csv(ERROR)}
    for N in (4, 11):
        base = next(r for r in mapped if int(r["N"]) == N
                    and float(r["salt_mM"]) == 1 and r["kind"] == "free")
        for variation in ("half_step", "double_L"):
            if (N, variation) in existing_error:
                continue
            row = calculate(spec, chemistry, N, 1, "free",
                            mesh_step_nm=float(base["mesh_step_nm"]) / 2 if variation == "half_step" else None,
                            outer_factor=2 if variation == "double_L" else 1)
            append_csv(ERROR, {"N": N, "variation": variation,
                               "sigma_base_C_m2": float(base["sigma_mean_C_m2"]),
                               "sigma_varied_C_m2": row["sigma_mean_C_m2"],
                               "absolute_change_C_m2": abs(row["sigma_mean_C_m2"] - float(base["sigma_mean_C_m2"])),
                               "status": row["status"]})
            print(f"N{N} {variation} {row['status']}", flush=True)
    errors = read_csv(ERROR)
    assert len(errors) == 4

    fixed_sources = {}
    for salt in (0.1, 1, 10):
        source = next(r for r in mapped if int(r["N"]) == 20
                      and float(r["salt_mM"]) == salt and r["kind"] == "free")
        fixed_sources[salt] = float(source["sigma_mean_C_m2"])
    existing_control = {(int(r["N"]), float(r["salt_mM"]), r["kind"]) for r in read_csv(CONTROL)}
    for salt in (0.1, 1, 10):
        for N in (1, 3, 11, 20):
            for kind in ("free", "backed"):
                if (N, float(salt), kind) in existing_control:
                    continue
                row = calculate(spec, chemistry, N, salt, kind, fixed_sigma=fixed_sources[salt])
                row.update({"fixed_sigma_source_C_m2": fixed_sources[salt],
                            "fixed_sigma_source": f"free N20 at c_s={salt:g} mM pH7; same reference for all N/kind",
                            "parameter_role": "MODEL_ONLY_FIXED_CHARGE_CONTROL"})
                append_csv(CONTROL, row)
                print(f"fixed N{N} c={salt:g} {kind} {row['status']}", flush=True)
    controlled = read_csv(CONTROL)
    assert len(controlled) == 24
    sigma4 = float(next(r["sigma_mean_C_m2"] for r in mapped if int(r["N"]) == 4
                        and float(r["salt_mM"]) == 1 and r["kind"] == "free"))
    sigma11 = float(next(r["sigma_mean_C_m2"] for r in mapped if int(r["N"]) == 11
                         and float(r["salt_mM"]) == 1 and r["kind"] == "free"))
    epsilon = {N: max(float(r["absolute_change_C_m2"]) for r in errors if int(r["N"]) == N)
               for N in (4, 11)}
    difference = abs(sigma4) - abs(sigma11)
    score = {"parameter_role": "MODEL_ONLY_LEGACY_REFERENCE_NOT_PE_FIT",
             "reduced_map_fields": len(mapped), "constant_charge_fields": len(controlled),
             "N4_N11_extra_fields": len(errors),
             "sigma_N4_C_m2": sigma4, "sigma_N11_C_m2": sigma11,
             "signed_sigma_N4_minus_N11_C_m2": sigma4 - sigma11,
             "magnitude_difference_C_m2": difference,
             "N4_error_C_m2": epsilon[4], "N11_error_C_m2": epsilon[11],
             "five_times_sum_error_C_m2": 5 * (epsilon[4] + epsilon[11]),
             "ten_percent_sigma_N11_C_m2": 0.1 * abs(sigma11),
             "model_only_enhancement_exceeds_both_thresholds": difference > max(0.1 * abs(sigma11),
                                                                               5 * (epsilon[4] + epsilon[11])),
             "not_material_validation": True}
    (OUT / "v9_MODEL_ONLY_精简图判定.json").write_text(json.dumps(score, ensure_ascii=False, indent=2) + "\n",
                                                        encoding="utf-8")
    print(json.dumps(score, ensure_ascii=False))


if __name__ == "__main__":
    main()
