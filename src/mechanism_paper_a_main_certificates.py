"""Closed-form checks for key C1/C2 W2 G1 E1 paper A branches."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0"
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
LP = ROOT / "data/processed/salt_transport/25_论文裁剪A五目标界_20261031.csv"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪A主条件解析核对_20261031.csv"
FIELDS = ("scenario", "window", "family", "delta_rho", "external",
          "B2600_over_B65_lower_closed", "R65_over_B65_lower_closed",
          "R2600_over_R65_lower_closed", "B2600_over_B65_upper_unbounded_witness",
          "R2600_over_B2600_lower_zero_witness",
          "R700_over_B700_upper_unbounded_witness",
          "max_relative_difference_from_LP", "scope")


def main() -> None:
    with (PACKAGE / "论文贡献可行性核查.json").open(encoding="utf-8-sig") as stream:
        length = json.load(stream)["length_bounds"]["G1"]
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    with LP.open(encoding="utf-8-sig", newline="") as stream:
        lp = list(csv.DictReader(stream))
    source_by_key = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    lp_by_key = {(r["scenario"], r["window"], r["target"]): r for r in lp
                 if r["family"] == "G1" and r["delta_rho"] == "0.2"
                 and r["external"] == "E1" and r["mobility_sensitivity_fraction"] == "0.0"}
    output = []
    for scenario in ("C1_T85", "C2_T85"):
        rho = {}
        for sample in ("B65", "B700", "B2600", "R65", "R700", "R2600"):
            row = source_by_key[(scenario, "W2", sample)]
            if row["status"] != "VALID" or float(row["lambda_s_inv"]) <= 0:
                raise ValueError("Invalid frozen rate")
            rho[sample] = 1 / float(row["lambda_s_inv"])
        lo = {s: .8 * value for s, value in rho.items()}
        hi = {s: 1.2 * value for s, value in rho.items()}
        thick = lo["B65"] / hi["B2600"] * length["2600"]["min_nm"] / length["65"]["max_nm"]
        treat = (lo["B65"] - hi["B2600"]) / hi["R65"] * length["65"]["min_nm"] / length["65"]["max_nm"]
        rthick = ((lo["R65"] - hi["R700"]) / (hi["R2600"] - hi["R700"])
                  * length["2600"]["min_nm"] / length["65"]["max_nm"])
        b_witness = hi["B2600"] < min(lo["B65"], lo["B700"])
        r_witness = hi["R700"] < min(lo["R65"], lo["R2600"])
        if not (b_witness and r_witness):
            raise ValueError("Unbounded branch construction failed")
        check = {"B2600_over_B65": thick, "R65_over_B65": treat,
                 "R2600_over_R65": rthick}
        errors = []
        for target, value in check.items():
            saved = float(lp_by_key[(scenario, "W2", target)]["lower_b_ratio"])
            errors.append(abs(value - saved) / value)
        if max(errors) > 1e-10:
            raise ValueError("Closed-form and LP lower bounds differ")
        output.append({
            "scenario": scenario, "window": "W2", "family": "G1",
            "delta_rho": .2, "external": "E1",
            "B2600_over_B65_lower_closed": thick,
            "R65_over_B65_lower_closed": treat,
            "R2600_over_R65_lower_closed": rthick,
            "B2600_over_B65_upper_unbounded_witness": b_witness,
            "R2600_over_B2600_lower_zero_witness": b_witness,
            "R700_over_B700_upper_unbounded_witness": r_witness,
            "max_relative_difference_from_LP": max(errors),
            "scope": "ARCHIVED_G1_OUTER_ENVELOPE_CONDITIONAL",
        })
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"main_conditions={len(output)} max_relative_difference="
          f"{max(r['max_relative_difference_from_LP'] for r in output):.3g}")


if __name__ == "__main__":
    main()
