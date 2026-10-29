"""Map each fixed core H-cell run to the planned quality-control evidence."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/experiment_runs.json"
OUT = ROOT / "data/processed/salt_transport/22_WP5核心运行质量证据关联_20261029.csv"
FIELDS = ("core_run_id", "group", "thickness_level", "parent_film_id", "PE_state",
          "salt_level", "C1_support_same_salt_run_ids", "C2_seal_same_salt_run_id",
          "C3_same_parent_state_blank_run_id", "C4_same_parent_state_salt_no_urea_run_id",
          "E3_same_parent_water_wait_run_id", "E4_same_parent_three_layer_run_id",
          "evidence_status")


def main() -> None:
    runs = json.loads(RUNS.read_text(encoding="utf-8-sig"))
    core = [r for r in runs if r["group"] in ("E1", "E2")]
    control = [r for r in runs if r["group"] not in ("E1", "E2")]
    if len(runs) != 54 or len(core) != 30 or len(control) != 24:
        raise ArithmeticError("Fixed H-cell matrix size changed")
    if any(r["status"] != "NOT_STARTED" for r in runs):
        raise ArithmeticError("Actual H-cell run status changed")
    output = []
    for run in core:
        same_salt = [c for c in control if c["salt_level"] == run["salt_level"]]
        c1 = sorted(c["run_id"] for c in same_salt if c["group"] == "C1")
        c2 = sorted(c["run_id"] for c in same_salt if c["group"] == "C2")
        c3 = sorted(c["run_id"] for c in control if c["group"] == "C3"
                    and c["parent_film_id"] == run["parent_film_id"]
                    and c["PE_state"] == run["PE_state"])
        c4 = sorted(c["run_id"] for c in control if c["group"] == "C4"
                    and c["parent_film_id"] == run["parent_film_id"]
                    and c["PE_state"] == run["PE_state"]
                    and c["salt_level"] == run["salt_level"])
        e3 = sorted(c["run_id"] for c in control if c["group"] == "E3"
                    and run["thickness_level"] == "thin" and run["salt_level"] == "high"
                    and c["parent_film_id"] == run["parent_film_id"])
        e4 = sorted(c["run_id"] for c in control if c["group"] == "E4"
                    and run["thickness_level"] == "thin" and run["salt_level"] == "high"
                    and run["PE_state"] == "E" and c["parent_film_id"] == run["parent_film_id"])
        if len(c1) != 3 or len(c2) != 1 or any(len(x) > 1 for x in (c3, c4, e3, e4)):
            raise ArithmeticError(f"Unexpected control coverage for {run['run_id']}")
        output.append({
            "core_run_id": run["run_id"], "group": run["group"],
            "thickness_level": run["thickness_level"],
            "parent_film_id": run["parent_film_id"], "PE_state": run["PE_state"],
            "salt_level": run["salt_level"],
            "C1_support_same_salt_run_ids": ";".join(c1),
            "C2_seal_same_salt_run_id": c2[0],
            "C3_same_parent_state_blank_run_id": ";".join(c3),
            "C4_same_parent_state_salt_no_urea_run_id": ";".join(c4),
            "E3_same_parent_water_wait_run_id": ";".join(e3),
            "E4_same_parent_three_layer_run_id": ";".join(e4),
            "evidence_status": "PLANNED_NOT_MEASURED",
        })
    if (sum(bool(r["C3_same_parent_state_blank_run_id"]) for r in output) != 12
            or sum(bool(r["C4_same_parent_state_salt_no_urea_run_id"]) for r in output) != 4
            or sum(bool(r["E3_same_parent_water_wait_run_id"]) for r in output) != 6
            or sum(bool(r["E4_same_parent_three_layer_run_id"]) for r in output) != 3):
        raise ArithmeticError("Unexpected limited-control coverage")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print("core=30 C1=30 C2=30 C3_same_parent=12 C4_P1=4 E3=6 E4=3")


if __name__ == "__main__":
    main()
