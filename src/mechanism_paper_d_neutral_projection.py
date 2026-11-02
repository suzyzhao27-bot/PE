"""Paper-scope D old-data-conditioned urea concentration envelopes, no ODE."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0"
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪D尿素浓度条件范围_20261031.csv"
GAIN_OUT = ROOT / "data/processed/salt_transport/25_论文裁剪D尿素处理响应范围_20261031.csv"
D_NA, D_CL, D_U, R_U_NM = 1.33e-9, 2.03e-9, 1.38e-9, .18
D_S = 2 * D_NA * D_CL / (D_NA + D_CL)
FIELDS = ("scenario", "sample_id", "h_nm", "time_min", "external_relation",
          "salt_exposure_min", "salt_exposure_max", "lambda_u_over_s_min",
          "lambda_u_over_s_max", "urea_receiver_min_mM", "urea_receiver_max_mM",
          "status", "scope")
GAIN_FIELDS = ("scenario", "h_nm", "window", "external_relation",
               "salt_gain_R_over_B", "urea_gain_min", "urea_gain_max",
               "status", "scope")


def hindrance(a_nm: int) -> float:
    x = R_U_NM / a_nm
    return (1 - x) ** 2 * (1 - 2.104 * x + 2.09 * x**3 - .95 * x**5)


def neutral_ratio_bounds(extrema: list[dict], h: int,
                         mobility_sensitivity_fraction: float = 0.) -> tuple[float, float]:
    ratios = []
    for e in extrema:
        a = e["a_nm"]
        if a not in (10, 25, 50):
            continue
        for tau in (1, 2, 4):
            if e["sigma_C_m2"] != 0 and tau * h / a < 5:
                continue
            entry = math.pi * a / (2 * tau)
            lu = h / hindrance(a) + entry
            ls_min = h / (e["Dbar_max_over_Ds"] * (1 + mobility_sensitivity_fraction)) + entry
            ls_max = h / (e["Dbar_min_over_Ds"] * (1 - mobility_sensitivity_fraction)) + entry
            ratios.append((D_U / D_S * min(1, ls_min / lu),
                           D_U / D_S * max(1, ls_max / lu)))
    return min(x[0] for x in ratios), max(x[1] for x in ratios)


def main() -> None:
    with (PACKAGE / "论文贡献可行性核查.json").open(encoding="utf-8-sig") as stream:
        pilot = json.load(stream)
    with (PACKAGE / "原矩阵论文用途子集.json").open(encoding="utf-8-sig") as stream:
        core = [x for x in json.load(stream) if x["paper_neutral_projection"]]
    if len(core) != 2 or {x["scenario"] for x in core} != {"C1_T85", "C2_T85"}:
        raise ValueError("Unexpected D scope")
    ratios = {h: neutral_ratio_bounds(pilot["operator_extrema"], h)
              for h in (65, 700, 2600)}
    for saved in pilot["analytic_neutral_ratio_diagnostic"]:
        actual = ratios[saved["h_nm"]]
        if not all(math.isclose(a, b, rel_tol=1e-12) for a, b in zip(
                actual, (saved["lambda_u_over_lambda_s_min"],
                         saved["lambda_u_over_lambda_s_max"]))):
            raise ValueError("Neutral ratio bound differs from archived pilot")

    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    if len(lookup) != len(source):
        raise ValueError("Duplicate source identity")
    output = []
    for case in core:
        scenario = case["scenario"]
        delta = case["delta_rho"]
        if delta != .2:
            raise ValueError("Unexpected D tolerance")
        for sample in ("B65", "B700", "B2600", "R65", "R700", "R2600"):
            h = int(sample[1:])
            qmin, qmax = ratios[h]
            for time_min, windows in ((60, (("W1", 3600),)),
                                      (100, (("W1", 3600), ("W2", 2400)))):
                zmin = zmax = 0.0
                for window, seconds in windows:
                    row = lookup[(scenario, window, sample)]
                    lam = float(row["lambda_s_inv"])
                    if row["status"] != "VALID" or lam <= 0:
                        raise ValueError("Invalid frozen window")
                    zmin += lam * seconds / (1 + delta)
                    zmax += lam * seconds / (1 - delta)
                common = dict(scenario=scenario, sample_id=sample, h_nm=h,
                              time_min=time_min, salt_exposure_min=zmin,
                              salt_exposure_max=zmax, scope="OLD_DEVICE_CONDITIONAL_TARGET")
                output.append(dict(common, external_relation="DIFFUSION_SCALED_SHARED",
                                   lambda_u_over_s_min=qmin, lambda_u_over_s_max=qmax,
                                   urea_receiver_min_mM=5 * -math.expm1(-qmin * zmin),
                                   urea_receiver_max_mM=5 * -math.expm1(-qmax * zmax),
                                   status="MODEL_ONLY_CONDITIONAL"))
                output.append(dict(common, external_relation="UREA_EXTERNAL_INDEPENDENT",
                                   lambda_u_over_s_min="NO_POSITIVE_BOUND",
                                   lambda_u_over_s_max="NO_FINITE_BOUND",
                                   urea_receiver_min_mM=0, urea_receiver_max_mM=5,
                                   status="ANALYTIC_CANCELLATION_CLOSURE"))
    if len(output) != 48:
        raise ValueError("Expected 48 D rows")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    gains = []
    for case in core:
        scenario = case["scenario"]
        delta = case["delta_rho"]
        for h in (65, 700, 2600):
            qmin, qmax = ratios[h]
            for window in ("W1", "W2"):
                b = lookup[(scenario, window, f"B{h}")]
                r = lookup[(scenario, window, f"R{h}")]
                if b["status"] != "VALID" or r["status"] != "VALID":
                    raise ValueError("Invalid gain input")
                salt_gain = float(r["lambda_s_inv"]) / float(b["lambda_s_inv"])
                common = dict(scenario=scenario, h_nm=h, window=window,
                              salt_gain_R_over_B=salt_gain,
                              scope="OLD_DEVICE_CONDITIONAL_GROUP_COMPARISON")
                gains.append(dict(common, external_relation="DIFFUSION_SCALED_SHARED",
                                  urea_gain_min=salt_gain * (qmin / qmax) * (1 - delta) / (1 + delta),
                                  urea_gain_max=salt_gain * (qmax / qmin) * (1 + delta) / (1 - delta),
                                  status="CONSERVATIVE_OUTER_ENVELOPE"))
                gains.append(dict(common, external_relation="UREA_EXTERNAL_INDEPENDENT",
                                  urea_gain_min="NO_POSITIVE_BOUND",
                                  urea_gain_max="NO_FINITE_BOUND",
                                  status="ANALYTIC_CANCELLATION"))
    if len(gains) != 24:
        raise ValueError("Expected 24 gain rows")
    with GAIN_OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=GAIN_FIELDS)
        writer.writeheader()
        writer.writerows(gains)
    print(f"concentration_rows={len(output)} gain_rows={len(gains)}")


if __name__ == "__main__":
    main()
