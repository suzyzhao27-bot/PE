"""Expand the fixed core H-cell matrix into mother-film B/E analysis pairs."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/experiment_runs.json"
OUT = ROOT / "data/processed/salt_transport/22_WP6母膜B_E配对运行入口_20261029.csv"
FIELDS = ("thickness_level", "salt_level", "salt_feed_plan_mM", "parent_block",
          "parent_film_id", "B_run_id", "E_run_id", "analysis_role",
          "B_measurement_status", "E_measurement_status", "note")


def main() -> None:
    runs = json.loads(RUNS.read_text(encoding="utf-8-sig"))
    if len(runs) != 54:
        raise ArithmeticError("Fixed H-cell matrix changed")
    core = [r for r in runs if r["group"] in ("E1", "E2")]
    if len(core) != 30:
        raise ArithmeticError("Expected 30 core runs")
    grouped = defaultdict(dict)
    for run in core:
        key = run["thickness_level"], run["salt_level"], run["replicate_block"]
        if run["PE_state"] in grouped[key]:
            raise ArithmeticError("Duplicate state in one mother-film block")
        grouped[key][run["PE_state"]] = run
    if len(grouped) != 15:
        raise ArithmeticError("Expected five conditions x three parent films")
    rows = []
    used_ids = []
    for (level, salt, block), pair in sorted(grouped.items()):
        if set(pair) != {"B", "E"}:
            raise ArithmeticError("Incomplete B/E parent pair")
        b, e = pair["B"], pair["E"]
        if b["parent_film_id"] != e["parent_film_id"] or b["salt_feed_mM"] != e["salt_feed_mM"]:
            raise ArithmeticError("Parent or salt mismatch")
        used_ids.extend((b["run_id"], e["run_id"]))
        rows.append({"thickness_level": level, "salt_level": salt,
                     "salt_feed_plan_mM": b["salt_feed_mM"], "parent_block": block,
                     "parent_film_id": b["parent_film_id"],
                     "B_run_id": b["run_id"], "E_run_id": e["run_id"],
                     "analysis_role": "CHARGE_CALIBRATION" if (level, salt) == ("thin", "low") else "HELD_OUT_SALT_PREDICTION",
                     "B_measurement_status": b["status"], "E_measurement_status": e["status"],
                     "note": "PAIR_BY_PARENT_NOT_BY_TIMEPOINT_OR_COUPON"})
    if len(used_ids) != 30 or len(set(used_ids)) != 30 or set(used_ids) != {r["run_id"] for r in core}:
        raise ArithmeticError("Core run coverage failed")
    if sum(r["analysis_role"] == "CHARGE_CALIBRATION" for r in rows) != 3:
        raise ArithmeticError("Calibration parent pair count changed")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"parent_pairs={len(rows)} core_runs={len(used_ids)} calibration_parent_pairs=3 "
          f"heldout_parent_pairs=12")


if __name__ == "__main__":
    main()
