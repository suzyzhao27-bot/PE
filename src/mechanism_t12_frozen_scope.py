"""Map ten unresolved PB targets to the actual frozen 54-pair examples."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = DATA / "22_T12_十项未决与封存示例范围_20261029.csv"
FIELDS = ("a_nm", "sigma_C_m2", "feed_mM", "relative_31_16",
          "frozen_charged_B_pair_count", "fixed_core_low_salt_candidate_pair_count",
          "scope_note")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    targets = [r for r in read("22_pb_target_recheck_20261028.csv")
               if r["status"] == "NUMERIC_UNRESOLVED_OVER_2_PERCENT"]
    curves = read("22_paired_illustrations_20261028.csv")
    if len(pairs) != 54 or len(targets) != 10 or len(curves) != 864:
        raise ArithmeticError("Unexpected fixed source size")
    charged = [p for p in pairs if p["sigma_B_C_m2"] != 0]
    if len(charged) != 36 or any(p["sigma_E_C_m2"] != 0 for p in pairs):
        raise ArithmeticError("Frozen charge-state scope changed")
    charged_curve_ids = {r["pair_id"] for r in curves if r["state"] == "B"
                         and r["probe"] == "salt" and r["numeric_status"] == "NUMERIC_UNRESOLVED_T12"}
    output = []
    for target in targets:
        a, sigma, feed = (float(target[k]) for k in ("a_nm", "sigma_C_m2", "feed_mM"))
        # The 54-pair generator fixes PoreTable(a_nm=25, sigma=-0.008).
        affected = [p for p in charged if a == 25 and p["sigma_B_C_m2"] == sigma
                    and p["c_feed_mM"] == feed]
        output.append({
            "a_nm": a, "sigma_C_m2": sigma, "feed_mM": feed,
            "relative_31_16": target["updated_31_16_relative_difference"],
            "frozen_charged_B_pair_count": len(affected),
            "fixed_core_low_salt_candidate_pair_count": sum(p["h_nm"] in (65, 2600)
                                                       for p in affected),
            "scope_note": ("CURRENT_FROZEN_54_PAIR_INPUT" if affected else
                           "NOT_USED_BY_CURRENT_54_PAIR_EXAMPLES_WP6_SCOPE_SEPARATE"),
        })
    affected_ids = {p["pair_id"] for p in charged if p["c_feed_mM"] == 10}
    if (sum(r["frozen_charged_B_pair_count"] for r in output) != 12
            or sum(r["fixed_core_low_salt_candidate_pair_count"] for r in output) != 8
            or charged_curve_ids != affected_ids):
        raise ArithmeticError("Frozen T12 label does not match target mapping")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print("unresolved_targets=10 used_by_current_54=1 charged_B_pairs=12 "
          "fixed_core_candidate_pairs=8")


if __name__ == "__main__":
    main()
