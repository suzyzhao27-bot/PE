"""Reproduce the frozen v8.0.1 unit and observation checks; no model fit."""

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/16_物理诊断与接续任务_v8.0.1"
if not PACKAGE.exists():
    PACKAGE = ROOT / "taskbook"  # standalone formal package copy
OUT = ROOT / "WP1"
MW = 58.44
C0 = 600.0  # mol/m3; declared protocol scenario, not measured for each film
DB = 2 * 1.33e-9 * 2.03e-9 / (1.33e-9 + 2.03e-9)
T_END = 120 * 60.0


def segment(sc):
    return "zero" if sc == 0 else "linear_0_10" if sc < 10 else "poly"


def concentration(sc):
    """Conditional NaCl concentration, mol/m3, from delta conductivity in uS/cm."""
    if sc < 0:
        raise ValueError("negative delta conductivity")
    if sc == 0:
        return 0.0
    at_10 = 10 ** (0.0165 + 0.9398 - 3.2793) * 1000 / MW
    if sc < 10:
        return sc * at_10 / 10
    x = math.log10(sc)
    return 10 ** (0.0165 * x * x + 0.9398 * x - 3.2793) * 1000 / MW


def conductivity(c):
    """Inverse of frozen conversion on its monotone branch, uS/cm."""
    if c < 0:
        raise ValueError("negative concentration")
    if c == 0:
        return 0.0
    join = concentration(10)
    if c < join:
        return 10 * c / join
    y = math.log10(c * MW / 1000)
    a, b, d = 0.0165, 0.9398, -3.2793 - y
    x = (-b + math.sqrt(b * b - 4 * a * d)) / (2 * a)
    return 10 ** x


def write_csv(path, rows, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def close(actual, expected, name):
    if not math.isclose(actual, expected, rel_tol=0.01, abs_tol=1e-14):
        raise AssertionError(f"{name}: {actual} != {expected}")


def main():
    salt = read_csv(PACKAGE / "inputs/salt/salt_transport_tidy.csv")
    roles = {r["curve_id"]: r["role"] for r in read_csv(PACKAGE / "inputs/salt/curve_roles_v7.csv")}
    assert len(salt) == 167 and len(roles) == 17
    groups = {}
    role_counts = {}
    for row in salt:
        assert row["conductivity_unit"] == "μS/cm"
        role = roles[row["curve_id"]]
        role_counts[role] = role_counts.get(role, 0) + 1
        groups.setdefault(row["curve_id"], []).append(row)
    assert role_counts == {"B": 42, "L": 70, "C": 14, "U": 20, "D": 21}
    endpoints = []
    for curve_id, points in groups.items():
        if roles[curve_id] != "B":
            continue
        assert len(points) == 7 and all(p["time_unit"] == "min" for p in points)
        points.sort(key=lambda p: float(p["time"]))
        h = int(points[0]["sample_id"].split()[0])
        assert float(points[0]["time"]) == 0 and float(points[-1]["time"]) == 120
        delta = float(points[-1]["conductivity"]) - float(points[0]["conductivity"])
        c = concentration(delta)
        fraction = c / C0
        lam = -math.log1p(-2 * fraction) / (2 * T_END)
        p_app = lam / 10.0  # alpha = 10 m^-1
        endpoints.append({
            "h_nm": h, "delta_uScm": delta, "conversion_segment": segment(delta),
            "c_eq_molm3": c, "receiver_fraction_C0": fraction,
            "lambda_app_s_inv": lam, "P_app_alpha10_m_s": p_app,
            "chi_app_alpha10": p_app * h * 1e-9 / DB,
            "h2_over_Db_s": (h * 1e-9) ** 2 / DB,
            "lag_upper_alpha10_epsilon035_s": 0.35 * h * 1e-9 / (4 * p_app),
            "time_evidence_level": "explicit_origin_units",
            "interpretation": "conditional endpoint-equivalent; not measured steady P",
        })
    endpoints.sort(key=lambda r: r["h_nm"])
    assert [r["h_nm"] for r in endpoints] == [65, 120, 200, 350, 700, 2600]
    audit = read_csv(PACKAGE / "audit/02_端点及必要量级.csv")
    audit_b = {int(r["h_nm"]): r for r in audit if r["role"] == "B"}
    for row in endpoints:
        ref = audit_b[row["h_nm"]]
        for own, theirs in (("delta_uScm", "delta_uScm"), ("c_eq_molm3", "c_eq_molm3"),
                            ("P_app_alpha10_m_s", "P_avg_at_alpha10_m_s"),
                            ("chi_app_alpha10", "chi_app"),
                            ("h2_over_Db_s", "tau_bulk_h2_D_s"),
                            ("lag_upper_alpha10_epsilon035_s", "lag_upper_if_Kle1_phi035_s")):
            close(row[own], float(ref[theirs]), f"{row['h_nm']} {own}")
    write_csv(OUT / "量纲与假设表.csv", endpoints, list(endpoints[0]))

    checks = []
    for sc in (0, 0.08, 1, 10, 100, 1363.25, 6872.54):
        c = concentration(sc)
        back = conductivity(c)
        err = abs(back - sc) / sc if sc else abs(back)
        assert err <= (1e-8 if sc else 1e-12)
        checks.append({"input_delta_uScm": sc, "conversion_segment": segment(sc),
                       "c_eq_molm3": c, "roundtrip_delta_uScm": back,
                       "relative_error_or_zero_absolute": err})
    write_csv(OUT / "七点电导往返记录.csv", checks, list(checks[0]))

    p = {r["h_nm"]: r for r in endpoints}
    chi_ratio = p[2600]["chi_app_alpha10"] / p[350]["chi_app_alpha10"]
    close(chi_ratio, 907.3752040917087, "2600/350 chi ratio")
    x_over_c = 1 / 0.01 - 0.01
    sigma = x_over_c * C0 * 96485.33212 * 25e-9 / 2
    close(sigma, 72.356762690091, "Donnan sigma")
    na, cl, f, gas, kelvin = 1.33e-9, 2.03e-9, 96485.33212, 8.314462618, 298.15
    ne_scale = f * f * (na + cl) / (gas * kelvin) * C0 * 1e4
    close(ne_scale, 75708.40927806833, "NE scale")
    molality = (35 / MW) / 0.965  # 3.5 wt% means 35 g salt per kg solution
    close(molality, 0.6206267976011889, "3.5 wt% molality")
    alpha_lag = {str(h): {str(t): 4 * p[h]["lambda_app_s_inv"] * t / (0.35 * h * 1e-9)
                          for t in (1200, 2400)} for h in (350, 2600)}
    supplement = json.loads((PACKAGE / "audit/06_审计补充量级.json").read_text(encoding="utf-8"))
    for h, values in alpha_lag.items():
        for t, value in values.items():
            close(value, supplement["lag_alpha_thresholds_phi035"][h][t], f"lag alpha {h} {t}")
    donor_remaining = 1 - p[2600]["receiver_fraction_C0"]
    gradient_remaining = 1 - 2 * p[2600]["receiver_fraction_C0"]
    assert sum(roles[r["curve_id"]] == "B" for r in salt) == 42
    summary = {
        "status": "v8_step1_existing_data_arithmetic_only",
        "source": "16_物理诊断与接续任务_v8.0.1/inputs/salt; frozen_spec_v8.json",
        "assumptions": "C0=600 mol/m3, two 10 mL baths, alpha=10 m^-1 endpoint example, 25 C conversion, epsilon=0.35 and K<=1 for lag bound",
        "Db_m2_s": DB, "NE_Q_uScm": ne_scale,
        "B2600_c_molm3": p[2600]["c_eq_molm3"],
        "B2600_receiver_fraction": p[2600]["receiver_fraction_C0"],
        "B2600_donor_remaining_fraction": donor_remaining,
        "B2600_concentration_gradient_remaining_fraction": gradient_remaining,
        "molality_3p5wt_mol_kg_water": molality,
        "alpha_1e5_effective_pool_depth_m": 1e-5,
        "alpha_lag_20_40min_epsilon035_m_inv": alpha_lag,
        "role_point_counts": role_counts,
        "alpha_1_10_100_area_cm2": [0.1, 1, 10],
        "lag_upper_alpha1_to100_all_B_s": [min(r["lag_upper_alpha10_epsilon035_s"] for r in endpoints) / 10,
                                               max(r["lag_upper_alpha10_epsilon035_s"] for r in endpoints) * 10],
        "chi_ratio_2600_350": chi_ratio,
        "sigma_C_m2_for_1pct_Donnan_a25nm": sigma,
        "conversion_join_c_molm3": concentration(10),
        "roundtrip_count": len(checks), "audit_B_endpoint_rows_compared": len(endpoints),
        "new_network_solver_run": False, "optimization_run": False,
    }
    (OUT / "v8_量级复核.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
