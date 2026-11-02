"""Audit links between frozen model comparisons and fixed H-cell runs."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/experiment_runs.json"
COMPARISONS = ROOT / "data/processed/salt_transport/22_W_C_WC_处理响应差距_20261028.csv"
OUT = ROOT / "data/processed/salt_transport/22_模型与实验条件接口_20261029.csv"
FIELDS = ("h_nm", "salt_feed_mM", "external", "thickness_level", "core_B_run_ids",
          "core_E_run_ids", "link_if_low_1mM", "link_if_low_10mM",
          "model_comparison_rows", "numeric_status", "interpretation")


def main() -> None:
    runs = json.loads(RUNS.read_text(encoding="utf-8-sig"))
    with COMPARISONS.open(encoding="utf-8-sig", newline="") as stream:
        comparisons = list(csv.DictReader(stream))
    if len(runs) != 54 or len(comparisons) != 162:
        raise ArithmeticError("Unexpected fixed matrix or comparison size")

    grouped = defaultdict(list)
    for row in comparisons:
        grouped[int(row["h_nm"]), int(row["salt_feed_mM"]), row["external"]].append(row)
    if len(grouped) != 18 or any(len(rows) != 9 for rows in grouped.values()):
        raise ArithmeticError("Incomplete model condition strata")

    level_by_h = {65: "thin", 700: "middle", 2600: "thick"}
    output = []
    for (h, feed, external), rows in sorted(grouped.items()):
        level = level_by_h[h]
        state_runs = {}
        for state in ("B", "E"):
            subset = [r for r in runs if r["group"] in ("E1", "E2")
                      and r["thickness_level"] == level and r["PE_state"] == state
                      and ((feed in (1, 10) and r["salt_level"] == "low")
                           or (feed == 600 and r["salt_level"] == "high"))]
            state_runs[state] = sorted(subset, key=lambda r: r["replicate_block"])
        linked = len(state_runs["B"]) == len(state_runs["E"]) == 3
        if linked and any(b["parent_film_id"] != e["parent_film_id"]
                          for b, e in zip(state_runs["B"], state_runs["E"])):
            raise ArithmeticError("Core B/E parents are not paired")
        if not linked and (state_runs["B"] or state_runs["E"]):
            raise ArithmeticError("Partial core link")
        if feed == 10 and any(r["numeric_status"] != "NUMERIC_UNRESOLVED_T12" for r in rows):
            raise ArithmeticError("Missing T12 warning")
        output.append({
            "h_nm": h, "salt_feed_mM": feed, "external": external,
            "thickness_level": level,
            "core_B_run_ids": ";".join(r["run_id"] for r in state_runs["B"]),
            "core_E_run_ids": ";".join(r["run_id"] for r in state_runs["E"]),
            "link_if_low_1mM": int(linked and feed != 10),
            "link_if_low_10mM": int(linked and feed != 1),
            "model_comparison_rows": len(rows),
            "numeric_status": "NUMERIC_UNRESOLVED_T12" if feed == 10 else "MODEL_ONLY_TARGET_CHECK",
            "interpretation": ("NO_FIXED_CORE_PAIR" if not linked else
                               "CONDITIONAL_NOMINAL_THICKNESS_EXTERNAL_NOT_MEASURED"),
        })
    for field in ("link_if_low_1mM", "link_if_low_10mM"):
        if sum(row[field] for row in output) != 10:
            raise ArithmeticError(f"Unexpected linked strata: {field}")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print("model_strata=18 comparisons=162 linked_per_low_choice=10 "
          "linked_comparisons_per_low_choice=90")


if __name__ == "__main__":
    main()
