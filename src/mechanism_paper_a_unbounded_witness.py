"""Construct denominator-zero limits for all fixed paper A upper bounds."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
BOUNDS = ROOT / "data/processed/salt_transport/25_论文裁剪A五目标界_20261031.csv"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪A上界无界构造_20261031.csv"
TARGETS = {"B2600_over_B65": ("B65", "B2600"),
           "R2600_over_R65": ("R65", "R2600"),
           "R65_over_B65": ("B65", "R65"),
           "R700_over_B700": ("B700", "R700"),
           "R2600_over_B2600": ("B2600", "R2600")}
GROUPS = {"B": ("B65", "B700", "B2600"), "R": ("R65", "R700", "R2600")}
FIELDS = ("case_id", "window", "target", "scenario", "Gamma", "delta_rho",
          "denominator_sample", "other_group_min_scaled_rho_upper_s",
          "denominator_rho_lower_s", "strict_witness_available",
          "LP_upper_status", "status_agrees", "constructive_limit")


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    with BOUNDS.open(encoding="utf-8-sig", newline="") as stream:
        bounds = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    if len(lookup) != len(source) or len(bounds) != 80:
        raise ValueError("Source or A table identity mismatch")
    output = []
    for row in bounds:
        numerator, denominator = TARGETS[row["target"]]
        delta, gamma = float(row["delta_rho"]), float(row["Gamma"])
        rho = {}
        for sample in GROUPS[denominator[0]]:
            frozen = lookup[(row["scenario"], row["window"], sample)]
            if frozen["status"] != "VALID" or float(frozen["lambda_s_inv"]) <= 0:
                raise ValueError("Invalid frozen rate")
            rho[sample] = 1 / float(frozen["lambda_s_inv"])
        other = [sample for sample in GROUPS[denominator[0]] if sample != denominator]
        r_den = (1 - delta) * rho[denominator]
        other_min = min((1 + delta) * rho[sample] * gamma for sample in other)
        # Set xi_den=1 and rho_den=rho_den^-. For other group mates,
        # set xi=Gamma and rho=rho^+. Then r_group tends to r_den from below.
        # If all other group totals exceed r_den, v_den -> 0 while v_num > 0.
        witness = other_min > r_den
        lp_unbounded = row["upper_status"] == "UNBOUNDED"
        output.append({
            "case_id": row["case_id"], "window": row["window"],
            "target": row["target"], "scenario": row["scenario"],
            "Gamma": gamma, "delta_rho": delta,
            "denominator_sample": denominator,
            "other_group_min_scaled_rho_upper_s": other_min,
            "denominator_rho_lower_s": r_den,
            "strict_witness_available": witness,
            "LP_upper_status": row["upper_status"],
            "status_agrees": witness == lp_unbounded,
            "constructive_limit": ("r_group -> denominator_rho_lower; v_denominator -> 0+"
                                   if witness else "NO_DENOMINATOR_ZERO_WITNESS_IN_THIS_FORM"),
        })
    mismatches = [r for r in output if not r["status_agrees"]]
    if mismatches:
        raise ArithmeticError(f"Upper-status construction mismatch: {len(mismatches)}")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(output)} constructed_unbounded="
          f"{sum(r['strict_witness_available'] for r in output)} mismatches=0")


if __name__ == "__main__":
    main()
