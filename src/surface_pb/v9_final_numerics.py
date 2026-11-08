"""Frozen 12-state final mesh/domain checks for fallback chemistry."""

import csv
import json
from pathlib import Path

from v9_geometry import family_geometry, pure_geometry
from v9_observation import solve_velocity
from v9_pb_1d import Chemical
from v9_pb_2d import solve_equilibrium

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP5/v9_MODEL_ONLY_十二状态数值复核.csv"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def save(rows):
    if not rows:
        return
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    decision = json.loads((ROOT / "WP4/v9_第4步判定.json").read_text(encoding="utf-8"))
    if decision["compatible_gamma_count"] != 0:
        raise RuntimeError("fallback chemistry only when no compatible calibration")
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    old = spec["chemistry"]["legacy_fallback"]
    chemistry = Chemical(old["Gamma_nm2"], old["pI"], old["b"])
    observations = {r["case_id"]: r for r in read_csv(ROOT / "WP1/v9_九组观测与条件盐度.csv")}
    technical = {r["case_id"]: r for r in read_csv(ROOT / "WP3/v9_四状态网格外域.csv")}
    map_rows = read_csv(ROOT / "WP5/v9_MODEL_ONLY_精简厚度图.csv")
    cases = []
    for draw in (182, 357, 604):
        for group in (3, 7, 10):
            source = observations[f"{draw}_pH{group}"]
            cases.append((f"{draw}_pH{group}_backed", family_geometry(spec, draw, "backed"),
                          float(source["conditional_NaCl_mM"]), float(source["pH_measured_mean"]),
                          float(source["H_full_gap_um_mean_assumed"])))
    for N in (1, 6, 20):
        cases.append((f"pure_N{N}_pH7_free", pure_geometry(spec, N, "free"), 1.0, 7.0, None))
    existing = {r["case_id"]: r for r in read_csv(OUT)} if OUT.exists() else {}
    extra_fields = 0
    for case_id, geometry, salt, pH, H in cases:
        if case_id in existing:
            continue
        reuse = technical.get(case_id)
        if reuse:
            sigma_base = float(reuse["sigma_base_C_m2"])
            sigma_refined = float(reuse["sigma_refined_C_m2"])
            sigma_outer = float(reuse["sigma_outer_C_m2"])
            zeta_base = float(reuse["zeta_base_mV"]) if H else None
            zeta_refined = float(reuse["zeta_refined_mV"]) if H else None
            zeta_outer = float(reuse["zeta_outer_mV"]) if H else None
            max_gauss = float(reuse["max_Gauss_rel"])
            max_species = float(reuse["max_species_M"])
            source = "REUSED_WP3_TECHNICAL_3_FIELDS"
        else:
            map_base = next((r for r in map_rows if case_id == "pure_N20_pH7_free"
                             and int(r["N"]) == 20 and float(r["salt_mM"]) == 1
                             and r["kind"] == "free"), None)
            base = None if map_base else solve_equilibrium(geometry, salt, pH, chemical=chemistry)
            extra_fields += 0 if map_base else 1
            step = float(map_base["mesh_step_nm"]) if map_base else base.mesh.target_nm
            refined = solve_equilibrium(geometry, salt, pH, chemical=chemistry,
                                        mesh_step_nm=step / 2)
            outer = solve_equilibrium(geometry, salt, pH, chemical=chemistry, outer_factor=2)
            extra_fields += 2
            sigma_base = float(map_base["sigma_mean_C_m2"]) if map_base else base.sigma_mean_C_m2
            sigma_refined = refined.sigma_mean_C_m2
            sigma_outer = outer.sigma_mean_C_m2
            if H:
                zeta_base = solve_velocity(geometry, base, H).zeta_app_mV
                zeta_refined = solve_velocity(geometry, refined, H).zeta_app_mV
                zeta_outer = solve_velocity(geometry, outer, H).zeta_app_mV
            else:
                zeta_base = zeta_refined = zeta_outer = None
            max_gauss = max([float(map_base["Gauss_relative"]) if map_base else base.Gauss_relative,
                             refined.Gauss_relative, outer.Gauss_relative])
            max_species = max([float(map_base["max_species_M"]) if map_base else base.max_species_M,
                               refined.max_species_M, outer.max_species_M])
            source = "REUSED_WP5_MAP_BASE_PLUS_2" if map_base else "NEW_BASE_REFINED_DOUBLE_L"
        epsilon_mesh = abs(sigma_refined - sigma_base)
        epsilon_outer = abs(sigma_outer - sigma_base)
        epsilon = max(epsilon_mesh, epsilon_outer)
        near_zero = abs(sigma_base) < 1e-5
        sigma_pass = (epsilon <= 1e-5 if near_zero else epsilon / abs(sigma_base) <= 0.02)
        zeta_mesh = abs(zeta_refined - zeta_base) if H else None
        zeta_outer_diff = abs(zeta_outer - zeta_base) if H else None
        zeta_pass = not H or max(zeta_mesh, zeta_outer_diff) <= 1
        physical = max_gauss <= 1e-6 and max_species <= 0.5
        row = {"case_id": case_id, "source": source,
               "parameter_role": "MODEL_ONLY_LEGACY_REFERENCE_NOT_PE_FIT",
               "salt_mM": salt, "pH": pH, "H_um": H if H else "NA",
               "sigma_base_C_m2": sigma_base, "sigma_refined_C_m2": sigma_refined,
               "sigma_double_L_C_m2": sigma_outer,
               "epsilon_mesh_C_m2": epsilon_mesh, "epsilon_double_L_C_m2": epsilon_outer,
               "epsilon_sigma_C_m2": epsilon,
               "epsilon_sigma_relative": epsilon / abs(sigma_base) if not near_zero else "NA",
               "zeta_base_mV": zeta_base if H else "NA",
               "zeta_refined_mV": zeta_refined if H else "NA",
               "zeta_double_L_mV": zeta_outer if H else "NA",
               "epsilon_zeta_mesh_mV": zeta_mesh if H else "NA",
               "epsilon_zeta_double_L_mV": zeta_outer_diff if H else "NA",
               "max_Gauss_relative": max_gauss, "max_species_M": max_species,
               "sigma_pass": sigma_pass, "zeta_pass": zeta_pass, "physical_pass": physical,
               "status": "PASS" if sigma_pass and zeta_pass and physical else "FAIL_RETAINED"}
        existing[case_id] = row
        save([existing[name] for name, *_ in cases if name in existing])
        print(f"{case_id}: {row['status']} ({source})", flush=True)
    rows = [existing[name] for name, *_ in cases]
    total_new_fields = sum(3 if r["source"] == "NEW_BASE_REFINED_DOUBLE_L" else
                           2 if r["source"] == "REUSED_WP5_MAP_BASE_PLUS_2" else 0 for r in rows)
    result = {"overall": "PASS" if all(r["status"] == "PASS" for r in rows) else "NOT_PASSED",
              "states": len(rows), "new_PB_fields_total": total_new_fields,
              "max_epsilon_sigma_relative": max(float(r["epsilon_sigma_relative"]) for r in rows
                                                if r["epsilon_sigma_relative"] != "NA"),
              "max_epsilon_zeta_mV": max(max(float(r["epsilon_zeta_mesh_mV"]),
                                          float(r["epsilon_zeta_double_L_mV"]))
                                          for r in rows if r["H_um"] != "NA"),
              "max_Gauss_relative": max(float(r["max_Gauss_relative"]) for r in rows),
              "max_species_M": max(float(r["max_species_M"]) for r in rows),
              "full_step5_status": "IN_PROGRESS_SENSITIVITIES_NOT_COMPLETE"}
    (ROOT / "WP5/v9_MODEL_ONLY_十二状态数值复核.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    if result["overall"] != "PASS":
        raise AssertionError("some frozen numerical states failed; inspect detailed CSV")


if __name__ == "__main__":
    main()
