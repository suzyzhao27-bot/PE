"""Focused validation of the 2D finite-PE sensitivity solver."""

import csv
import json
from pathlib import Path

from v9_geometry import family_geometry, pure_geometry
from v9_observation import solve_velocity
from v9_pb_1d import Chemical, F, R, T
from v9_pb_2d import solve_equilibrium
from v9_pb_dielectric_2d import solve_dielectric

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    old = spec["chemistry"]["legacy_fallback"]
    chemistry = Chemical(old["Gamma_nm2"], old["pI"], old["b"])
    N1 = pure_geometry(spec, 1, "free")
    zero = solve_dielectric(N1, 1, 7, Chemical(0, 4, -1), mesh_step_nm=3)
    zero_psi = max(max(abs(zero.u)), max(abs(zero.solid_u))) * R * T / F * 1000
    main = solve_equilibrium(N1, 1, 7, chemical=chemistry, mesh_step_nm=3)
    small_pe = solve_dielectric(N1, 1, 7, chemistry, eps_PE_r=0.001,
                                mesh_step_nm=3)
    gauss_limit_sigma_rel = (abs(small_pe.sigma_mean_C_m2 - main.sigma_mean_C_m2)
                             / abs(main.sigma_mean_C_m2))
    with (ROOT / "WP1/v9_九组观测与条件盐度.csv").open(encoding="utf-8-sig", newline="") as stream:
        observation = next(r for r in csv.DictReader(stream) if r["case_id"] == "182_pH7")
    geometry = family_geometry(spec, 182, "backed")
    salt = float(observation["conditional_NaCl_mM"])
    pH = float(observation["pH_measured_mean"])
    H = float(observation["H_full_gap_um_mean_assumed"])
    base = solve_dielectric(geometry, salt, pH, chemistry)
    refined = solve_dielectric(geometry, salt, pH, chemistry,
                               mesh_step_nm=base.mesh.target_nm / 2)
    outer = solve_dielectric(geometry, salt, pH, chemistry, outer_factor=2)
    zeta = [solve_velocity(geometry, field, H).zeta_app_mV for field in (base, refined, outer)]
    sigma_change = max(abs(refined.sigma_mean_C_m2 - base.sigma_mean_C_m2),
                       abs(outer.sigma_mean_C_m2 - base.sigma_mean_C_m2))
    sigma_change_rel = sigma_change / abs(base.sigma_mean_C_m2)
    zeta_change = max(abs(zeta[1] - zeta[0]), abs(zeta[2] - zeta[0]))
    result = {"zero_psi_abs_mV": zero_psi,
              "zero_sigma_C_m2": zero.sigma_mean_C_m2,
              "small_epsilon_PE_r": 0.001,
              "small_epsilon_limit_sigma_relative_error": gauss_limit_sigma_rel,
              "representative_case": "182_pH7_backed",
              "epsilon_PE_r": 2.3,
              "sigma_base_C_m2": base.sigma_mean_C_m2,
              "sigma_refined_C_m2": refined.sigma_mean_C_m2,
              "sigma_double_L_C_m2": outer.sigma_mean_C_m2,
              "max_sigma_change_relative": sigma_change_rel,
              "zeta_base_mV": zeta[0], "zeta_refined_mV": zeta[1],
              "zeta_double_L_mV": zeta[2],
              "max_zeta_change_mV": zeta_change,
              "max_Gauss_relative": max(x.Gauss_relative for x in (base, refined, outer)),
              "max_species_M": max(x.max_species_M for x in (base, refined, outer)),
              "additional_refinement_PB_fields": 2,
              "status": "PASS" if zero_psi <= 1e-6 and zero.sigma_mean_C_m2 == 0
                                  and gauss_limit_sigma_rel <= 1e-4
                                  and sigma_change_rel <= 0.02 and zeta_change <= 1
                                  and all(x.Gauss_relative <= 1e-6 and not x.physical_PB_warning
                                          for x in (base, refined, outer)) else "NOT_PASSED"}
    (ROOT / "WP5/v9_二维介电求解器核对.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    if result["status"] != "PASS":
        raise AssertionError("finite dielectric solver validation did not pass")


if __name__ == "__main__":
    main()
