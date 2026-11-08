"""Lock the prescribed two calibration and eight held-out core conditions to H runs."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/experiment_runs.json"
OUT = ROOT / "data/processed/salt_transport/22_WP6校准与留出实验行_20261029.csv"
FIELDS = ("thickness_level", "salt_level", "salt_feed_plan_mM", "PE_state",
          "source_group", "run_ids_P1_P2_P3", "parent_film_ids_P1_P2_P3",
          "analysis_role", "measured_status")


def main() -> None:
    runs = json.loads(RUNS.read_text(encoding="utf-8-sig"))
    if len(runs) != 54 or any(r["status"] != "NOT_STARTED" for r in runs):
        raise ArithmeticError("Fixed run matrix changed or experiments started")
    core = [r for r in runs if r["group"] in ("E1", "E2")]
    if len(core) != 30:
        raise ArithmeticError("Expected 30 core H-cell runs")
    grouped = defaultdict(list)
    for run in core:
        grouped[run["thickness_level"], run["salt_level"], run["PE_state"]].append(run)
    if len(grouped) != 10:
        raise ArithmeticError("Expected ten core condition means")
    output = []
    for (level, salt, state), subset in sorted(grouped.items()):
        subset.sort(key=lambda r: r["replicate_block"])
        if [r["replicate_block"] for r in subset] != [1, 2, 3]:
            raise ArithmeticError("Missing independent P1/P2/P3 mother-film block")
        if len({r["group"] for r in subset}) != 1:
            raise ArithmeticError("Mixed source groups")
        role = "CHARGE_CALIBRATION_WITH_INDEPENDENT_UREA_B" if level == "thin" and salt == "low" else "HELD_OUT_SALT_PREDICTION"
        output.append({
            "thickness_level": level, "salt_level": salt,
            "salt_feed_plan_mM": subset[0]["salt_feed_mM"],
            "PE_state": state, "source_group": subset[0]["group"],
            "run_ids_P1_P2_P3": ";".join(r["run_id"] for r in subset),
            "parent_film_ids_P1_P2_P3": ";".join(r["parent_film_id"] for r in subset),
            "analysis_role": role, "measured_status": "NOT_STARTED",
        })
    if sum(r["analysis_role"].startswith("CHARGE_CALIBRATION") for r in output) != 2:
        raise ArithmeticError("Calibration split is not two conditions")
    if sum(r["analysis_role"] == "HELD_OUT_SALT_PREDICTION" for r in output) != 8:
        raise ArithmeticError("Held-out split is not eight conditions")
    by_condition = {(r["thickness_level"], r["salt_level"], r["PE_state"]): r for r in output}
    for level, salt, _ in grouped:
        if (level, salt, "B") not in by_condition or (level, salt, "E") not in by_condition:
            raise ArithmeticError("Missing B/E condition pair")
        b = by_condition[level, salt, "B"]["parent_film_ids_P1_P2_P3"]
        e = by_condition[level, salt, "E"]["parent_film_ids_P1_P2_P3"]
        if b != e:
            raise ArithmeticError("B/E mother-film mismatch")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print("core_runs=30 condition_means=10 calibration=2 held_out=8 mother_film_pairs=5x3")


if __name__ == "__main__":
    main()
