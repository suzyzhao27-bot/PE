"""Formal 1D PB reference checks toward, but not completing, v9 gate 3."""

import csv
import json
import math
from pathlib import Path

from fibre_pb import Adsorption, solve as old_fibre_solve
from v9_analytic import state
from v9_pb_1d import Chemical, EPS0, F, R, T, solve

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP3"


def main():
    data = json.loads((ROOT / "WP1/v9_数据审计.json").read_text(encoding="utf-8"))
    algebra = json.loads((ROOT / "WP2/v9_解析核对.json").read_text(encoding="utf-8"))
    assert data["status"] == "PASS_DATA_AUDIT" and algebra["status"] == "PASS_ANALYTIC_ONLY"
    zero = solve(fixed_sigma=0.0)
    zero_chemical = solve(chemical=Chemical(0.0, 4.0, -1.0))
    zero_abs_mV = max(abs(float(zero.surface_psi_mV)),
                      abs(float(zero_chemical.surface_psi_mV)))
    assert zero_abs_mV <= 1e-6 and zero.gauss_absolute == 0
    positive, negative = solve(fixed_sigma=0.005), solve(fixed_sigma=-0.005)
    sigma_scale = math.sqrt(8 * EPS0 * 78.3 * R * T * 1.0)
    psi_grahame_mV = 2 * R * T / F * math.asinh(0.005 / sigma_scale) * 1000
    grahame_rel = max(abs(float(positive.surface_psi_mV) - psi_grahame_mV),
                      abs(float(negative.surface_psi_mV) + psi_grahame_mV)) / psi_grahame_mV
    sign_abs = abs(float(positive.surface_psi_mV + negative.surface_psi_mV))
    chemical_reversal = abs(state(0.7, 4.0, 0.05, 3.0, -0.5)[3]
                            + state(-0.7, 4.0, 0.05, 5.0, -0.5)[3])
    assert grahame_rel <= 0.01 and sign_abs <= 1e-8 and chemical_reversal <= 1e-14
    legacy = Chemical(0.05, 4.372113695091497, -0.8682652885468274)
    new = solve(kind="radial", inner_nm=13.3, outer_nm=39.5, n=512,
                salt_mM=1.0, pH=7.0, chemical=legacy, outer_condition="cell")
    old = old_fibre_solve(salt_mM=1.0, pH=7.0, radius_nm=13.3,
                          cell_nm=39.5, n=512,
                          adsorption=Adsorption(0.05, {"H": -20, "OH": -50}))
    sigma_rel = abs(new.sigma_C_m2 - old["surface_sigma_C_m2"]) / abs(old["surface_sigma_C_m2"])
    psi_rel = abs(float(new.surface_psi_mV) - old["surface_psi_mV"]) / abs(old["surface_psi_mV"])
    assert sigma_rel <= 0.01 and psi_rel <= 0.01
    assert max(positive.gauss_relative, negative.gauss_relative, new.gauss_relative) <= 1e-6
    assert min(zero.min_species_M, positive.min_species_M,
               negative.min_species_M, new.min_species_M) >= 0

    rows = [
        (1, "PASS_DATA_ONLY", "9 groups, 36 ramps", "WP1/v9_数据审计.json"),
        (2, "PASS_ANALYTIC_ONLY", f"max dσ/du rel {algebra['max_derivative_relative_error']:.3g}", "WP2/v9_解析核对.json"),
        (3, "PARTIAL_PB_ZERO", f"planar ψ={zero_abs_mV:.3g} mV; ζ not implemented", "this run"),
        (4, "PASS_1D_REFERENCE", f"Grahame rel {grahame_rel:.6g}; sign abs {sign_abs:.3g} mV", "this run"),
        (5, "PASS_1D_REFERENCE", f"v3 σ rel {sigma_rel:.3g}; ψ rel {psi_rel:.3g}", "this run"),
        (6, "NOT_RUN", "finite dielectric slab implementation pending", ""),
        (7, "PARTIAL_1D_ONLY", f"1D Gauss max {max(positive.gauss_relative,negative.gauss_relative,new.gauss_relative):.3g}; finite fiber pending", "this run"),
        (8, "NOT_RUN", "velocity and ζ observation pending", ""),
        (9, "NOT_RUN", "four finite-fiber mesh states pending", ""),
        (10, "NOT_RUN", "four finite-fiber outer-domain states pending", ""),
    ]
    with (OUT / "v9_十项技术进度.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["check_id", "status", "observed_metric_or_gap", "evidence"])
        writer.writerows(rows)
    record = {"overall": "TECHNICAL_NOT_PASSED", "PB_1D_reference_run": True,
              "finite_fiber_PB_run": False, "calibration_run": False,
              "zero_potential_abs_mV": zero_abs_mV,
              "Grahame_surface_potential_rel_error": grahame_rel,
              "fixed_charge_sign_symmetry_abs_mV": sign_abs,
              "chemistry_state_reversal_abs_C_m2": chemical_reversal,
              "v3_radial_sigma_rel_error": sigma_rel,
              "v3_radial_psi_rel_error": psi_rel,
              "max_1D_Gauss_relative": max(positive.gauss_relative, negative.gauss_relative, new.gauss_relative),
              "one_dimensional_results_do_not_validate_finite_fiber_or_zeta": True}
    (OUT / "v9_一维参考核对.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, ensure_ascii=False))


if __name__ == "__main__":
    main()
