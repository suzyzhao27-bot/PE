"""Compare frozen charged-salt curves with a 30-second RK4 step, without new PB calls."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from mechanism_transport import PoreTable, run_curve


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
FROZEN = DATA / "22_paired_illustrations_20261028.csv"
OUT = DATA / "22_带电盐时步收敛_60对30秒_20261029.csv"
FIELDS = ("pair_id", "h_nm", "salt_feed_mM", "mechanism", "external", "state",
          "largest_receiver_difference_uM", "largest_relative_difference_to_30s",
          "time_of_largest_absolute_difference_min", "max_mass_relative_error_30s",
          "numeric_status")


def main() -> None:
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    with FROZEN.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(pairs) != 54 or len(rows) != 864:
        raise ArithmeticError("Unexpected fixed source size")
    frozen = {(r["pair_id"], r["state"], r["probe"], int(float(r["time_min"]))): r for r in rows}
    if len(frozen) != 864:
        raise ArithmeticError("Duplicate frozen state")
    table = PoreTable()
    output = []
    for pair in pairs:
        sigma = float(pair["sigma_B_C_m2"])
        if sigma == 0:
            continue
        initial = float(pair["c_feed_mM"])
        curve = run_curve(initial, "salt", float(pair["h_nm"]), float(pair["b_B"]),
                          sigma, pair["external"], table, step_seconds=30.0)
        gaps = []
        relative = []
        for item in curve:
            minute = int(item["time_min"])
            old = frozen[pair["pair_id"], "B", "salt", minute]
            if abs(item["donor_mol_m3"] + item["receiver_mol_m3"] - initial) > 1e-8 * initial:
                raise ArithmeticError("30-second reservoir concentration sum changed")
            if minute:
                gap = abs(item["receiver_mol_m3"] - float(old["receiver_mol_m3"]))
                gaps.append((minute, gap))
                relative.append(gap / item["receiver_mol_m3"] if item["receiver_mol_m3"] else 0.0)
        minute, gap = max(gaps, key=lambda x: x[1])
        output.append({"pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                       "salt_feed_mM": pair["c_feed_mM"], "mechanism": pair["mechanism"],
                       "external": pair["external"], "state": "B",
                       "largest_receiver_difference_uM": gap * 1000,
                       "largest_relative_difference_to_30s": max(relative),
                       "time_of_largest_absolute_difference_min": minute,
                       "max_mass_relative_error_30s": max(r["mass_relative_error"] for r in curve),
                       "numeric_status": "TIMESTEP_COMPARISON_MODEL_ONLY"})
        print(f'{pair["pair_id"]} {pair["c_feed_mM"]}mM', flush=True)
    if len(output) != 36:
        raise ArithmeticError("Expected 36 charged-salt curves")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"curves={len(output)} max_difference_uM={max(r['largest_receiver_difference_uM'] for r in output):.8g} "
          f"max_relative={max(r['largest_relative_difference_to_30s'] for r in output):.8g} "
          f"over_2pct={sum(r['largest_relative_difference_to_30s'] > .02 for r in output)}")


if __name__ == "__main__":
    main()
