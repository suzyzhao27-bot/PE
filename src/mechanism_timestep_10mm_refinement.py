"""Refine the twelve charged 10 mM curves from 30 to 15 seconds."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from mechanism_transport import PoreTable, run_curve


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = DATA / "22_10mM带电盐时步细化_30对15秒_20261029.csv"
FIELDS = ("pair_id", "h_nm", "mechanism", "external", "time_min", "receiver_30s_uM",
          "receiver_15s_uM", "difference_uM", "relative_difference_to_15s",
          "mass_relative_error_15s", "numeric_status")


def main() -> None:
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    if len(pairs) != 54:
        raise ArithmeticError("Unexpected fixed matrix")
    table = PoreTable()
    output = []
    selected = 0
    for pair in pairs:
        if float(pair["c_feed_mM"]) != 10 or float(pair["sigma_B_C_m2"]) == 0:
            continue
        selected += 1
        args = (10.0, "salt", float(pair["h_nm"]), float(pair["b_B"]),
                float(pair["sigma_B_C_m2"]), pair["external"], table)
        coarse = run_curve(*args, step_seconds=30.0)
        fine = run_curve(*args, step_seconds=15.0)
        for a, b in zip(coarse, fine, strict=True):
            minute = int(a["time_min"])
            if minute == 0:
                continue
            if abs(b["donor_mol_m3"] + b["receiver_mol_m3"] - 10) > 1e-7:
                raise ArithmeticError("Reservoir concentration sum changed")
            gap = abs(a["receiver_mol_m3"] - b["receiver_mol_m3"])
            output.append({"pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                           "mechanism": pair["mechanism"], "external": pair["external"],
                           "time_min": minute, "receiver_30s_uM": a["receiver_mol_m3"] * 1000,
                           "receiver_15s_uM": b["receiver_mol_m3"] * 1000,
                           "difference_uM": gap * 1000,
                           "relative_difference_to_15s": gap / b["receiver_mol_m3"],
                           "mass_relative_error_15s": b["mass_relative_error"],
                           "numeric_status": "TIMESTEP_COMPARISON_MODEL_ONLY"})
        print(pair["pair_id"], flush=True)
    if selected != 12 or len(output) != 36:
        raise ArithmeticError("Expected 12 pairs and 36 nonzero observations")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"pairs={selected} rows={len(output)} max_difference_uM={max(r['difference_uM'] for r in output):.8g} "
          f"max_relative={max(r['relative_difference_to_15s'] for r in output):.8g}")


if __name__ == "__main__":
    main()
