"""Frozen v9 charge-regulation algebra and scale checks; no PB field solve."""

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP2"
E = 1.602176634e-19
F = 96485.33212
R = 8.314462618
T = 298.15


def state(u, pH, gamma_nm2, pI, b):
    """Return occupancies, sigma, d sigma/du, and positive regulation capacitance."""
    a = math.log(10) * (b + pI - pH) - u
    c = math.log(10) * (b - pI + pH) + u
    anchor = max(0.0, a, c)
    w0, wp, wm = math.exp(-anchor), math.exp(a - anchor), math.exp(c - anchor)
    denominator = w0 + wp + wm
    theta0, thetap, thetam = (w0 / denominator, wp / denominator, wm / denominator)
    mean_z = thetap - thetam
    variance = thetap + thetam - mean_z * mean_z
    charge_capacity = E * gamma_nm2 * 1e18
    sigma = charge_capacity * mean_z
    derivative = -charge_capacity * variance
    c_reg = -derivative * F / (R * T)
    return theta0, thetap, thetam, sigma, derivative, c_reg, variance


def save(name, records):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    fallback = spec["chemistry"]["legacy_fallback"]
    gamma_grid = spec["chemistry"]["Gamma_nm2_grid"]
    assert len(gamma_grid) == 7 and spec["physics"]["T_K"] == T
    chemistry = []
    for gamma in gamma_grid:
        for pH in (3.0, 7.0, 10.0):
            for u in (-2.0, 0.0, 2.0):
                t0, tp, tm, sigma, analytic, c_reg, variance = state(
                    u, pH, gamma, fallback["pI"], fallback["b"])
                # At nearly saturated occupancy a 1e-6 step loses digits to
                # subtraction; 1e-4 stays below the frozen 1e-5 error gate.
                step = 1e-4
                central = (state(u + step, pH, gamma, fallback["pI"], fallback["b"])[3]
                           - state(u - step, pH, gamma, fallback["pI"], fallback["b"])[3]) / (2 * step)
                relative = abs(central - analytic) / max(abs(analytic), 1e-15)
                assert abs(t0 + tp + tm - 1) < 1e-14
                assert abs(sigma) <= E * gamma * 1e18 * (1 + 1e-14)
                assert analytic <= 0 and c_reg >= 0 and relative < 1e-5
                chemistry.append({"Gamma_nm2": gamma, "pH": pH, "u": u,
                                  "theta0": t0, "theta_plus": tp, "theta_minus": tm,
                                  "sigma_C_m2": sigma, "capacity_abs_C_m2": E * gamma * 1e18,
                                  "variance_z": variance, "d_sigma_du_C_m2": analytic,
                                  "finite_difference_d_sigma_du_C_m2": central,
                                  "derivative_rel_error": relative,
                                  "C_reg_F_m2": c_reg,
                                  "parameter_status": "LEGACY_ILLUSTRATION_NOT_FITTED"})
    save("v9_化学占据与导数核对.csv", chemistry)

    with (ROOT / "WP1/v9_九组观测与条件盐度.csv").open(encoding="utf-8-sig", newline="") as stream:
        observed = list(csv.DictReader(stream))
    ld = [float(row["conditional_lambda_D_nm"]) for row in observed]
    g = [float(spec["geometry"]["family"][str(draw)]["g_nm"]) for draw in (182, 357, 604)]
    small_field_min = 1 / math.tanh(max(g) / (2 * min(ld))) - 1
    small_field_max = 1 / math.tanh(min(g) / (2 * max(ld))) - 1
    pure = spec["geometry"]["pure_h"]
    d_nm, g_nm = pure["d_nm"], pure["g_nm"]
    scales = []
    for n in (1, 4, 11, 20, 40):
        h_nm = d_nm + (n - 1) * (d_nm + g_nm)
        scales.append({"N": n, "h_pure_nm": h_nm, "d_fixed_nm": d_nm,
                       "g_fixed_nm": g_nm, "edge_layer_fraction_2_over_N": 2 / n if n > 2 else "NA",
                       "h_over_min_conditional_lambda_D": h_nm / min(ld),
                       "h_over_max_conditional_lambda_D": h_nm / max(ld),
                       "meaning": "geometry/algebra only; no PB charge prediction"})
    save("v9_纯厚度与屏蔽尺度.csv", scales)
    check = {"status": "PASS_ANALYTIC_ONLY", "PB_solver_run": False,
             "charge_regulation_cases": len(chemistry),
             "max_derivative_relative_error": max(row["derivative_rel_error"] for row in chemistry),
             "lambda_D_conditional_nm_range": [min(ld), max(ld)],
             "dry_projected_gap_nm_range": [min(g), max(g)],
             "linear_parallel_wall_relative_potential_correction_range": [small_field_min, small_field_max],
             "pure_N4_h_nm": scales[1]["h_pure_nm"], "pure_N11_h_nm": scales[2]["h_pure_nm"],
             "symmetric_dense_PE_slab_single_face_sigma_h_dependence": "zero under identical baths and chemistry, symmetric stable solution",
             "parameter_status": "legacy fallback used for algebraic illustration only"}
    (OUT / "v9_解析核对.json").write_text(json.dumps(check, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(check, ensure_ascii=False))


if __name__ == "__main__":
    main()
