"""Audit all 24 fixed H-cell controls and their valid core-run comparisons."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/experiment_runs.json"
OUT = ROOT / "data/processed/salt_transport/22_WP5控制运行与核心关联_20261029.csv"
FIELDS = ("control_run_id", "group", "control_state", "salt_level",
          "parent_film_id", "linked_core_run_ids", "link_basis",
          "replicate_role", "measured_status")


def main() -> None:
    runs = json.loads(RUNS.read_text(encoding="utf-8-sig"))
    if len(runs) != 54 or {r["run_id"] for r in runs} != {f"H{i:03d}" for i in range(1, 55)}:
        raise ArithmeticError("Fixed 54-run identity changed")
    if any(r["status"] != "NOT_STARTED" for r in runs):
        raise ArithmeticError("Actual run status changed")
    core = [r for r in runs if r["group"] in ("E1", "E2")]
    controls = [r for r in runs if r["group"] not in ("E1", "E2")]
    if len(core) != 30 or len(controls) != 24:
        raise ArithmeticError("Unexpected core/control split")
    output = []
    for control in controls:
        group = control["group"]
        linked = []
        if group == "E3":
            linked = [r for r in core if r["thickness_level"] == "thin"
                      and r["salt_level"] == "high"
                      and r["parent_film_id"] == control["parent_film_id"]]
            basis, role = "SAME_PARENT_THIN_HIGH_B_AND_E", "WATER_WAIT_CONTROL_NOT_MODEL_W"
        elif group == "E4":
            linked = [r for r in core if r["thickness_level"] == "thin"
                      and r["salt_level"] == "high" and r["PE_state"] == "E"
                      and r["parent_film_id"] == control["parent_film_id"]]
            basis, role = "SAME_PARENT_SINGLE_LAYER_E_HIGH", "THREE_LAYER_SERIAL_CONTROL"
        elif group == "C3":
            linked = [r for r in core if r["thickness_level"] == "thin"
                      and r["PE_state"] == control["PE_state"]
                      and r["parent_film_id"] == control["parent_film_id"]]
            basis, role = "SAME_PARENT_AND_STATE_LOW_PLUS_HIGH", "ZERO_SOLUTE_BACKGROUND_NOT_MATERIAL_REPLICATE"
        elif group == "C4":
            linked = [r for r in core if r["thickness_level"] == "thin"
                      and r["PE_state"] == control["PE_state"]
                      and r["salt_level"] == control["salt_level"]
                      and r["parent_film_id"] == control["parent_film_id"]]
            basis, role = "SAME_PARENT_STATE_SALT_WITH_UREA", "UREA_PERTURBATION_SCREEN_P1_ONLY"
        elif group == "C1":
            basis, role = "SALT_LEVEL_SUPPORT_REFERENCE_NO_ONE_TO_ONE_FILM", "EXTERNAL_RESISTANCE_REFERENCE_NOT_MATERIAL_REPLICATE"
        elif group == "C2":
            basis, role = "SALT_LEVEL_SEAL_REFERENCE_NO_ONE_TO_ONE_FILM", "LEAKAGE_BACKGROUND_NOT_MATERIAL_REPLICATE"
        else:
            raise ArithmeticError(f"Unexpected control group {group}")
        expected = {"E3": 2, "E4": 1, "C1": 0, "C2": 0, "C3": 2, "C4": 1}[group]
        if len(linked) != expected:
            raise ArithmeticError(f"Unexpected control link count: {control['run_id']}")
        output.append({
            "control_run_id": control["run_id"], "group": group,
            "control_state": control["PE_state"], "salt_level": control["salt_level"],
            "parent_film_id": control["parent_film_id"],
            "linked_core_run_ids": ";".join(sorted(r["run_id"] for r in linked)),
            "link_basis": basis, "replicate_role": role,
            "measured_status": "NOT_STARTED",
        })
    if {r["control_run_id"] for r in output} != {f"H{i:03d}" for i in range(31, 55)}:
        raise ArithmeticError("Control run ID coverage mismatch")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print("controls=24 E3=3 E4=3 C1=6 C2=2 C3=6 C4=4 exact_links=16")


if __name__ == "__main__":
    main()
