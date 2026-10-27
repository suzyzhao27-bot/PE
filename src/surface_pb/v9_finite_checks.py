"""Frozen v9 finite-fiber checks, retaining failures and unrun dielectric gate."""

import csv
import json
import math
import time
from pathlib import Path

from v9_geometry import family_geometry, pure_geometry
from v9_observation import empty_slit_half_reference, flat_reference_zeta, solve_velocity
from v9_pb_1d import Chemical, F, R, T, Z, solve as solve_1d
from v9_pb_2d import solve_equilibrium

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP3"


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    old = spec["chemistry"]["legacy_fallback"]
    chemistry = Chemical(old["Gamma_nm2"], old["pI"], old["b"])
    with (ROOT / "WP1/v9_九组观测与条件盐度.csv").open(encoding="utf-8-sig", newline="") as stream:
        observations = {r["case_id"]: r for r in csv.DictReader(stream)}
    reference = json.loads((OUT / "v9_一维参考核对.json").read_text(encoding="utf-8"))
    assert reference["overall"] == "TECHNICAL_NOT_PASSED"
    checks = []
    cases = [("182_pH7_backed", family_geometry(spec, 182, "backed"), observations["182_pH7"]),
             ("604_pH3_backed", family_geometry(spec, 604, "backed"), observations["604_pH3"]),
             ("604_pH10_backed", family_geometry(spec, 604, "backed"), observations["604_pH10"]),
             ("pure_N1_pH7_free", pure_geometry(spec, 1, "free"), None)]
    total_started = time.perf_counter()
    for name, geometry, obs in cases:
        salt = float(obs["conditional_NaCl_mM"]) if obs else 1.0
        pH = float(obs["pH_measured_mean"]) if obs else 7.0
        H = float(obs["H_full_gap_um_mean_assumed"]) if obs else None
        record = {"case_id": name, "kind": geometry.kind, "c_s_mM": salt, "pH": pH,
                  "H_um": H if H else "NA", "parameter_status": "LEGACY_CODE_REFERENCE_NOT_FITTED"}
        started = time.perf_counter()
        try:
            base = solve_equilibrium(geometry, salt, pH, chemical=chemistry)
            refined = solve_equilibrium(geometry, salt, pH, chemical=chemistry,
                                        mesh_step_nm=base.mesh.target_nm / 2)
            outer = solve_equilibrium(geometry, salt, pH, chemical=chemistry,
                                      outer_factor=2)
            zeta_base = zeta_refined = zeta_outer = None
            flow_residual = top_spread = None
            if H:
                ob = solve_velocity(geometry, base, H)
                of = solve_velocity(geometry, refined, H)
                oo = solve_velocity(geometry, outer, H)
                zeta_base, zeta_refined, zeta_outer = ob.zeta_app_mV, of.zeta_app_mV, oo.zeta_app_mV
                flow_residual = max(ob.flow_residual_scaled, of.flow_residual_scaled, oo.flow_residual_scaled)
                top_spread = max(ob.v_top_relative_spread, of.v_top_relative_spread, oo.v_top_relative_spread)
            sigma_mesh_abs = abs(refined.sigma_mean_C_m2 - base.sigma_mean_C_m2)
            sigma_outer_abs = abs(outer.sigma_mean_C_m2 - base.sigma_mean_C_m2)
            sigma_mesh_rel = sigma_mesh_abs / max(abs(base.sigma_mean_C_m2), 1e-30)
            sigma_outer_rel = sigma_outer_abs / max(abs(base.sigma_mean_C_m2), 1e-30)
            zeta_mesh_abs = abs(zeta_refined - zeta_base) if H else None
            zeta_outer_abs = abs(zeta_outer - zeta_base) if H else None
            near_zero = abs(base.sigma_mean_C_m2) < 1e-5
            sigma_gate = (max(sigma_mesh_abs, sigma_outer_abs) <= 1e-5 if near_zero
                          else max(sigma_mesh_rel, sigma_outer_rel) <= 0.02)
            zeta_gate = (not H or max(zeta_mesh_abs, zeta_outer_abs) <= 1.0)
            physics_gate = all(x.Gauss_relative is not None and x.Gauss_relative <= 1e-6
                               and x.min_species_M >= 0 and not x.physical_PB_warning
                               for x in (base, refined, outer))
            top_u = max(abs(base.u[base.mesh.outer_cell]))
            top_mobile = max(abs(Z @ base.concentrations_mol_m3[:, base.mesh.outer_cell]))
            record.update({"base_cells": len(base.u), "refined_cells": len(refined.u),
                           "outer_cells": len(outer.u), "base_step_nm": base.mesh.target_nm,
                           "base_L_nm": base.mesh.L_nm,
                           "sigma_base_C_m2": base.sigma_mean_C_m2,
                           "sigma_refined_C_m2": refined.sigma_mean_C_m2,
                           "sigma_outer_C_m2": outer.sigma_mean_C_m2,
                           "sigma_mesh_abs_C_m2": sigma_mesh_abs,
                           "sigma_outer_abs_C_m2": sigma_outer_abs,
                           "sigma_mesh_rel": sigma_mesh_rel,
                           "sigma_outer_rel": sigma_outer_rel,
                           "zeta_base_mV": zeta_base if H else "NA",
                           "zeta_refined_mV": zeta_refined if H else "NA",
                           "zeta_outer_mV": zeta_outer if H else "NA",
                           "zeta_mesh_abs_mV": zeta_mesh_abs if H else "NA",
                           "zeta_outer_abs_mV": zeta_outer_abs if H else "NA",
                           "max_Gauss_rel": max(x.Gauss_relative for x in (base, refined, outer)),
                           "max_species_M": max(x.max_species_M for x in (base, refined, outer)),
                           "max_flow_residual_scaled": flow_residual if H else "NA",
                           "max_v_top_relative_spread": top_spread if H else "NA",
                           "base_top_psi_abs_mV": top_u * (R * T / F) * 1000,
                           "base_top_net_ion_mol_m3": top_mobile,
                           "sigma_gate": sigma_gate, "zeta_gate": zeta_gate,
                           "physics_gate": physics_gate,
                           "status": "PASS_NON_DIELECTRIC_SUBSET" if sigma_gate and zeta_gate and physics_gate else "FAIL_RETAINED"})
        except Exception as exc:
            record.update({"status": "NUMERICALLY_UNRESOLVED", "error": f"{type(exc).__name__}: {exc}"})
        record["elapsed_s"] = time.perf_counter() - started
        checks.append(record)
        print(f"{name}: {record['status']} ({record['elapsed_s']:.2f} s)", flush=True)
    columns = list(dict.fromkeys(k for row in checks for k in row))
    with (OUT / "v9_四状态网格外域.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(checks)
    empty_rel, empty_min = empty_slit_half_reference()
    flat = solve_1d(fixed_sigma=-0.005)
    hs = flat_reference_zeta(flat, 100)
    hs_rel = abs(hs - flat.surface_psi_mV) / abs(flat.surface_psi_mV)
    assert empty_rel <= 0.01 and empty_min >= 0 and hs_rel <= 0.01
    all_pass = all(row["status"] == "PASS_NON_DIELECTRIC_SUBSET" for row in checks)
    technical_rows = [
        (1, "PASS_DATA_ONLY", "9 groups and 36 retained ramps"),
        (2, "PASS_ANALYTIC_ONLY", "63 occupancy/derivative examples"),
        (3, "PARTIAL_ZERO", "1D potential zero; 2D zero PB passed previously; integrated ζ zero requires formal run"),
        (4, "PASS_1D_REFERENCE", f"Grahame rel {reference['Grahame_surface_potential_rel_error']:.6g}"),
        (5, "PASS_1D_REFERENCE", f"v3 σ/ψ rel {reference['v3_radial_sigma_rel_error']:.3g}/{reference['v3_radial_psi_rel_error']:.3g}"),
        (6, "NOT_RUN", "finite epsilon_PE=2.3 coupled dielectric slab pending"),
        (7, "PASS_NON_DIELECTRIC_SUBSET" if all_pass else "FAIL_RETAINED", "finite free/backed Gauss/site/nonnegativity on four states"),
        (8, "PASS_OBSERVATION_REFERENCE", f"empty Poiseuille rel {empty_rel:.3g}; planar HS rel {hs_rel:.3g}"),
        (9, "PASS_FOUR_STATES" if all_pass else "FAIL_RETAINED", "four base/refined comparisons in detailed table"),
        (10, "PASS_FOUR_STATES" if all_pass else "FAIL_RETAINED", "four base/double-L comparisons in detailed table"),
    ]
    with (OUT / "v9_十项技术进度.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["check_id", "status", "observed_metric_or_gap"])
        writer.writerows(technical_rows)
    summary = {"overall": "TECHNICAL_NOT_PASSED",
               "finite_fiber_four_states_all_pass_non_dielectric_subset": all_pass,
               "finite_dielectric_check6_run": False,
               "empty_slit_Poiseuille_rel_error": empty_rel,
               "planar_HS_rel_error": hs_rel,
               "states": {r["case_id"]: r["status"] for r in checks},
               "total_elapsed_s": time.perf_counter() - total_started,
               "calibration_run": False}
    (OUT / "v9_二维技术进度.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
