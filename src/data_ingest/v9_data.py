"""Audit the frozen v9 electrokinetic inputs without fitting or solving PB."""

import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
INPUT = PACKAGE / "inputs/zeta"
OUT = ROOT / "WP1"
EPS0 = 8.8541878128e-12
EPS_W = 78.3 * EPS0
F = 96485.33212
R = 8.314462618
T = 298.15
DIFFUSION = {"Na": 1.33e-9, "Cl": 2.03e-9, "H": 9.31e-9, "OH": 5.27e-9}


def rows(name):
    with (INPUT / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, records):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def close(actual, expected, label, tolerance=1e-9):
    if not math.isclose(float(actual), float(expected), rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{label}: {actual} != {expected}")


def conditional_bath(pH, conductivity_mS_m):
    lam = {name: F * F * value / (R * T) for name, value in DIFFUSION.items()}
    c_h = 1000 * 10 ** (-pH)
    c_oh = 1000 * 10 ** (pH - 14)
    extra_na = max(c_oh - c_h, 0)
    extra_cl = max(c_h - c_oh, 0)
    c_s = ((conductivity_mS_m / 1000) - lam["H"] * c_h
           - lam["OH"] * c_oh - lam["Na"] * extra_na
           - lam["Cl"] * extra_cl) / (lam["Na"] + lam["Cl"])
    if c_s <= 0:
        raise AssertionError("nonpositive conditional background NaCl")
    ionic_strength = c_s + max(c_h, c_oh)
    lambda_d_nm = 1e9 * math.sqrt(EPS_W * R * T / (2 * F * F * ionic_strength))
    return c_s, ionic_strength, lambda_d_nm


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    assert spec["primary_data"]["groups"] == 9 and spec["primary_data"]["pressure_ramps"] == 36
    originals = json.loads((INPUT / "original_independent_cells.json").read_text(encoding="utf-8"))
    workbooks = {int(w["workbook_name"].removesuffix(".xls")): w["sheets"][0]["values"] for w in originals}
    ramps = rows("zeta_ramps36.csv")
    summaries = rows("zeta_summary.csv")
    baths = {(int(r["draw"]), int(r["pH_group"])): r for r in rows("conditional_bath9.csv")}
    assert len(workbooks) == 3 and len(ramps) == 36 and len(summaries) == len(baths) == 9
    assert set(workbooks) == {182, 357, 604}

    checked_ramps = []
    group_rows = []
    max_hs_taskbook_delta = 0.0
    max_hs_report_delta = 0.0
    max_summary_delta = 0.0
    ramp_keys = set()
    for summary in summaries:
        draw, label = int(summary["draw_ratio_label"]), int(summary["pH_nominal_label"])
        key = (draw, label)
        sheet = workbooks[draw]
        summary_row = int(summary["summary_row_1based"])
        original_summary = sheet[summary_row - 1]
        matching = sorted((r for r in ramps if (int(r["draw"]), int(r["pH_group"])) == key),
                          key=lambda r: int(r["ramp"]))
        assert len(matching) == 4 and [int(r["ramp"]) for r in matching] == [1, 2, 3, 4]
        for name, col in (("pH_mean", 0), ("conductivity_mean_mS_per_m", 1),
                          ("zeta_mean_mV_reported", 3), ("zeta_stdev_mV_reported_ddof0", 4)):
            close(summary[name], original_summary[col], f"{key} summary {name}")
        for r in matching:
            row_num = int(r["source_row_1based"])
            assert row_num == summary_row + int(r["ramp"]) + 1
            original = sheet[row_num - 1]
            for name, col in (("ramp", 0), ("dI_dp_nA_mbar", 1), ("I_R2", 2),
                              ("dQ_dp_ml_min_mbar_signed", 3), ("Q_R2", 4),
                              ("pH", 5), ("conductivity_mS_m", 6), ("gap_um", 7),
                              ("zeta_reported_mV", 8)):
                close(r[name], original[col], f"{key} ramp {r['ramp']} {name}")
            ramp_key = (draw, label, int(r["ramp"]))
            assert ramp_key not in ramp_keys
            ramp_keys.add(ramp_key)
            dI_over_dQ = float(r["dI_dp_nA_mbar"]) * 0.06 / abs(float(r["dQ_dp_ml_min_mbar_signed"]))
            # The delivered reconstruction uses the reported instrument-current
            # orientation. Section 02 writes a leading minus sign but does not
            # certify that orientation; retain this discrepancy in the audit.
            hs = (float(r["gap_um"]) * 1e-6) ** 2 * dI_over_dQ / (12 * EPS_W) * 1000
            delta = hs - float(r["zeta_reported_mV"])
            max_hs_taskbook_delta = max(max_hs_taskbook_delta, abs(hs - float(r["zeta_HS_reconstructed_mV"])))
            max_hs_report_delta = max(max_hs_report_delta, abs(delta))
            checked_ramps.append({"draw": draw, "pH_group": label, "ramp": r["ramp"],
                                  "source_row_1based": row_num, "dI_dp_nA_mbar": r["dI_dp_nA_mbar"],
                                  "I_R2": r["I_R2"],
                                  "dQ_dp_ml_min_mbar_signed": r["dQ_dp_ml_min_mbar_signed"],
                                  "Q_R2": r["Q_R2"], "pH": r["pH"],
                                  "conductivity_mS_m": r["conductivity_mS_m"],
                                  "gap_um": r["gap_um"], "zeta_reported_mV": r["zeta_reported_mV"],
                                  "zeta_HS_mV_audit": hs, "HS_minus_reported_mV": delta,
                                  "flow_reverse": float(r["dQ_dp_ml_min_mbar_signed"]) < 0,
                                  "retained_even_if_low_R2": True})
        values = {name: [float(r[name]) for r in matching] for name in
                  ("pH", "conductivity_mS_m", "gap_um", "zeta_reported_mV")}
        means = {name: statistics.mean(seq) for name, seq in values.items()}
        sd = statistics.stdev(values["zeta_reported_mV"])
        for name, calculated in (("pH_mean", means["pH"]),
                                  ("conductivity_mean_mS_per_m", means["conductivity_mS_m"]),
                                  ("zeta_mean_mV_recomputed", means["zeta_reported_mV"]),
                                  ("zeta_sd_four_ramps_ddof1_mV", sd)):
            old = float(summary[name])
            max_summary_delta = max(max_summary_delta, abs(calculated - old))
            close(calculated, old, f"{key} {name}")
        bath = baths[key]
        cs, ionic, ld = conditional_bath(means["pH"], means["conductivity_mS_m"])
        for field, calculated in (("conditional_NaCl_mM", cs),
                                  ("conditional_I_mM", ionic), ("lambda_D_nm", ld)):
            close(calculated, bath[field], f"{key} {field}", tolerance=1e-6)
        close(means["gap_um"], bath["gap_um_mean"], f"{key} H")
        close(means["zeta_reported_mV"], bath["zeta_mV"], f"{key} zeta")
        group_rows.append({"case_id": f"{draw}_pH{label}", "draw": draw,
                           "pH_group": label, "pH_measured_mean": means["pH"],
                           "conductivity_mS_m_mean": means["conductivity_mS_m"],
                           "H_full_gap_um_mean_assumed": means["gap_um"],
                           "zeta_reported_mean_mV": means["zeta_reported_mV"],
                           "zeta_ramp_sd_ddof1_mV": sd,
                           "comparison_scale_mV": max(5.0, sd),
                           "conditional_NaCl_mM": cs, "conditional_I_mM": ionic,
                           "conditional_lambda_D_nm": ld,
                           "role": "calibration" if draw == 182 else "retrospective_prediction",
                           "independent_film_n": "unknown", "ramp_n": 4})
    assert len(ramp_keys) == 36 and len(group_rows) == 9
    assert all(r["independent_membrane_n"] == "unknown" for r in summaries)
    assert max_hs_taskbook_delta < 1e-7
    assert max_hs_report_delta < 5.0
    assert max_summary_delta < 1e-9
    assert Counter(int(r["draw"]) for r in ramps) == {182: 12, 357: 12, 604: 12}
    OUT.mkdir(exist_ok=True)
    write_csv("v9_36条扫描审计.csv", checked_ramps)
    write_csv("v9_九组观测与条件盐度.csv", group_rows)
    result = {"status": "PASS_DATA_AUDIT", "model_run": False, "fit_run": False,
              "workbooks": 3, "groups": 9, "ramps": 36,
              "reverse_flow_ramps": sum(r["flow_reverse"] for r in checked_ramps),
              "min_I_R2": min(float(r["I_R2"]) for r in checked_ramps),
              "min_Q_R2": min(float(r["Q_R2"]) for r in checked_ramps),
              "max_abs_HS_minus_reported_mV": max_hs_report_delta,
              "max_abs_HS_minus_taskbook_mV": max_hs_taskbook_delta,
              "max_abs_summary_recompute_difference": max_summary_delta,
              "conditional_NaCl_range_mM": [min(r["conditional_NaCl_mM"] for r in group_rows),
                                              max(r["conditional_NaCl_mM"] for r in group_rows)],
              "independent_film_n": "unknown", "H_full_gap": "assumed",
              "HS_sign_convention": "matches delivered reconstruction with reported current sign and abs(flow slope); section 02 displayed formula has opposite leading sign; instrument orientation unresolved",
              "bath_identity": "conditional NaCl+HCl/NaOH, not measured composition"}
    (OUT / "v9_数据审计.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
