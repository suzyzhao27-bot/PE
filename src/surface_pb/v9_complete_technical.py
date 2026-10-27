"""Close the ten frozen v9 technical checks; no calibration or material claim."""

import csv
import json
import math
from pathlib import Path

from v9_dielectric_slab import solve_slab
from v9_geometry import family_geometry
from v9_observation import empty_slit_half_reference, flat_reference_zeta, solve_velocity
from v9_pb_1d import Chemical, F, R, T, solve as solve_1d
from v9_pb_2d import solve_equilibrium

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP3"
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    data = json.loads((ROOT / "WP1/v9_数据审计.json").read_text(encoding="utf-8"))
    algebra = json.loads((ROOT / "WP2/v9_解析核对.json").read_text(encoding="utf-8"))
    reference = json.loads((OUT / "v9_一维参考核对.json").read_text(encoding="utf-8"))
    finite = read_csv(OUT / "v9_四状态网格外域.csv")
    obs = {r["case_id"]: r for r in read_csv(ROOT / "WP1/v9_九组观测与条件盐度.csv")}

    zero_case = obs["182_pH7"]
    g = family_geometry(spec, 182, "backed")
    zero_pb = solve_equilibrium(g, float(zero_case["conditional_NaCl_mM"]),
                                float(zero_case["pH_measured_mean"]), fixed_sigma=0.0)
    zero_ob = solve_velocity(g, zero_pb, float(zero_case["H_full_gap_um_mean_assumed"]))
    zero_psi_mV = float(max(abs(zero_pb.u)) * R * T / F * 1000)
    zero_pass = (zero_psi_mV <= 1e-6 and zero_pb.sigma_mean_C_m2 == 0
                 and abs(zero_ob.zeta_app_mV) <= 1e-6
                 and abs(zero_ob.I_half_A) <= 1e-20)

    old = spec["chemistry"]["legacy_fallback"]
    chemistry = Chemical(old["Gamma_nm2"], old["pI"], old["b"])
    slabs = [solve_slab(h, chemistry, chemistry) for h in (20, 200, 2000)]
    sigma = [s.sigma_left_C_m2 for s in slabs]
    sigma_range_rel = (max(sigma) - min(sigma)) / max(abs(sigma[0]), 1e-30)
    asym = solve_slab(20, Chemical(old["Gamma_nm2"], 4.2, old["b"]),
                      Chemical(old["Gamma_nm2"], 4.8, old["b"]))
    dielectric_pass = (sigma_range_rel <= 0.01
                       and all(abs(s.sigma_left_C_m2 - s.sigma_right_C_m2) <= 1e-12
                               and abs(s.E_PE_V_m) <= 1e-8
                               and s.total_Gauss_absolute_C_m2 <= 1e-6 for s in slabs)
                       and abs(asym.E_PE_V_m) > 1 and abs(asym.displacement_PE_C_m2) > 1e-8
                       and asym.total_Gauss_absolute_C_m2 <= 1e-6)

    empty_rel, empty_min = empty_slit_half_reference()
    positive = solve_1d(fixed_sigma=0.005)
    negative = solve_1d(fixed_sigma=-0.005)
    hs_pos = flat_reference_zeta(positive, 100)
    hs_neg = flat_reference_zeta(negative, 100)
    hs_rel = max(abs(hs_pos - positive.surface_psi_mV) / abs(positive.surface_psi_mV),
                 abs(hs_neg - negative.surface_psi_mV) / abs(negative.surface_psi_mV))
    observation_pass = (empty_rel <= 0.01 and empty_min >= 0 and hs_rel <= 0.01
                        and hs_pos > 0 > hs_neg and zero_pass)

    four = (len(finite) == 4 and all(r["status"] == "PASS_NON_DIELECTRIC_SUBSET"
                                      for r in finite))
    finite_physics = (four and all(float(r["max_Gauss_rel"]) <= 1e-6
                                   and float(r["max_species_M"]) <= 0.5
                                   and r["physics_gate"] == "True" for r in finite))
    mesh = four and all(r["sigma_gate"] == "True" and r["zeta_gate"] == "True"
                        for r in finite)
    outer = (mesh and all(float(r["base_top_psi_abs_mV"]) < 1e-3
                          and float(r["base_top_net_ion_mol_m3"]) < 1e-3
                          and (r["kind"] == "free" or float(r["max_v_top_relative_spread"]) < 1e-6)
                          for r in finite))
    # The upper bounds here are diagnostic checks of decay, not fitted model parameters.
    rows = [
        (1, data["status"] == "PASS_DATA_AUDIT" and data["groups"] == 9 and data["ramps"] == 36,
         f"{data['groups']} groups, {data['ramps']} ramps"),
        (2, algebra["status"] == "PASS_ANALYTIC_ONLY" and algebra["max_derivative_relative_error"] <= 1e-5,
         f"63 examples, derivative relative error {algebra['max_derivative_relative_error']:.6g}"),
        (3, zero_pass and reference["zero_potential_abs_mV"] <= 1e-6,
         f"2D zero ψ={zero_psi_mV:.6g} mV, σ={zero_pb.sigma_mean_C_m2:.6g}, ζ={zero_ob.zeta_app_mV:.6g} mV"),
        (4, reference["Grahame_surface_potential_rel_error"] <= 0.01
         and reference["fixed_charge_sign_symmetry_abs_mV"] <= 1e-8
         and reference["chemistry_state_reversal_abs_C_m2"] <= 1e-14,
         f"Grahame relative {reference['Grahame_surface_potential_rel_error']:.6g}, reversal checked"),
        (5, reference["v3_radial_sigma_rel_error"] <= 0.01 and reference["v3_radial_psi_rel_error"] <= 0.01,
         f"v3 σ/ψ relative {reference['v3_radial_sigma_rel_error']:.6g}/{reference['v3_radial_psi_rel_error']:.6g}"),
        (6, dielectric_pass,
         f"εPE=2.3 dense slab h=20/200/2000 nm; σ range relative {sigma_range_rel:.6g}; symmetric E_PE=0; asymmetric E_PE={asym.E_PE_V_m:.6g} V/m"),
        (7, finite_physics,
         f"four free/backed states; max Gauss relative {max(float(r['max_Gauss_rel']) for r in finite):.6g}; face capacity asserted in solver"),
        (8, observation_pass, f"Poiseuille relative {empty_rel:.6g}; ±HS relative {hs_rel:.6g}; zero current checked"),
        (9, mesh, f"four states; max Δσ relative {max(float(r['sigma_mesh_rel']) for r in finite):.6g}; max Δζ {max(float(r['zeta_mesh_abs_mV']) for r in finite if r['kind'] == 'backed'):.6g} mV"),
        (10, outer, f"four states L doubled; top |ψ| max {max(float(r['base_top_psi_abs_mV']) for r in finite):.6g} mV; top velocity spread max {max(float(r['max_v_top_relative_spread']) for r in finite if r['kind'] == 'backed'):.6g}"),
    ]
    for name in ("v9_十项技术验收.csv", "v9_十项技术进度.csv"):
        with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["check_id", "status", "observed_evidence"])
            writer.writerows((i, "PASS" if passed else "FAIL", evidence) for i, passed, evidence in rows)
    result = {
        "overall": "TECHNICAL_PASS" if all(passed for _, passed, _ in rows) else "TECHNICAL_NOT_PASSED",
        "passed_checks": sum(int(passed) for _, passed, _ in rows),
        "total_checks": 10,
        "calibration_run": False,
        "dense_slab_scope": "1D coupled finite dielectric only; 2D fiber dielectric sensitivity remains step 5",
        "dielectric": [s.__dict__ for s in slabs],
        "asymmetric_dielectric_probe": asym.__dict__,
        "zero_case": {"psi_abs_mV": zero_psi_mV, "sigma_C_m2": zero_pb.sigma_mean_C_m2,
                      "zeta_mV": zero_ob.zeta_app_mV, "I_half_A": zero_ob.I_half_A},
        "top_tail": {r["case_id"]: {"psi_abs_mV": float(r["base_top_psi_abs_mV"]),
                                     "net_ion_mol_m3": float(r["base_top_net_ion_mol_m3"])}
                     for r in finite},
    }
    (OUT / "v9_技术验收.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                                           encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("overall", "passed_checks", "total_checks")}, ensure_ascii=False))
    if result["overall"] != "TECHNICAL_PASS":
        raise AssertionError("v9 frozen technical gates not all passed; inspect CSV")


if __name__ == "__main__":
    main()
