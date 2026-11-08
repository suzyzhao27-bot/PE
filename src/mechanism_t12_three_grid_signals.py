"""Combine existing 16/31-node and new 35-node charged trajectory outputs."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_T12_10mM三网格接收信号范围_20261029.csv"
FIELDS = ("pair_id", "h_nm", "mechanism", "external", "time_min",
          "receiver_16_mM", "receiver_31_mM", "receiver_35_mM",
          "grid_spread_uM", "grid_spread_relative_to_31",
          "fixed_core_candidate_if_low_10mM", "numeric_status")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    old = read("22_T12_10mM完整时序差_20261029.csv")
    new = read("22_pb_postfreeze_10mM完整时序影响_20261029.csv")
    if len(old) != 48 or len(new) != 48:
        raise ArithmeticError("Unexpected curve size")
    refined = {(r["pair_id"], int(r["time_min"])): r for r in new}
    if len(refined) != 48:
        raise ArithmeticError("Duplicate post-freeze state")
    output = []
    for row in old:
        key = row["pair_id"], int(row["time_min"])
        other = refined[key]
        if any(row[k] != other[k] for k in ("mechanism", "external")) or float(row["h_nm"]) != float(other["h_nm"]):
            raise ArithmeticError("Curve identity mismatch")
        c16, c31, c35 = (float(row["receiver_16_mM"]), float(row["receiver_31_mM"]),
                         float(other["augmented_receiver_mM"]))
        if abs(c31 - float(other["frozen_receiver_mM"])) > 1e-10:
            raise ArithmeticError("Original 31-node value differs")
        spread = max(c16, c31, c35) - min(c16, c31, c35)
        output.append({
            "pair_id": key[0], "h_nm": row["h_nm"],
            "mechanism": row["mechanism"], "external": row["external"],
            "time_min": key[1], "receiver_16_mM": c16,
            "receiver_31_mM": c31, "receiver_35_mM": c35,
            "grid_spread_uM": spread * 1000,
            "grid_spread_relative_to_31": spread / c31 if c31 else 0.0,
            "fixed_core_candidate_if_low_10mM": int(float(row["h_nm"]) in (65, 2600)),
            "numeric_status": "INTERPOLATION_GRID_SPREAD_T12_UNRESOLVED",
        })
    if len(output) != 48 or sum(r["fixed_core_candidate_if_low_10mM"] for r in output) != 32:
        raise ArithmeticError("Unexpected candidate scope")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    later = [r for r in output if r["time_min"]]
    print(f"states=48 noninitial=36 max_spread_uM={max(r['grid_spread_uM'] for r in later):.8g} "
          f"max_relative={max(r['grid_spread_relative_to_31'] for r in later):.8g} "
          f"over_2pct={sum(r['grid_spread_relative_to_31'] > .02 for r in later)}")


if __name__ == "__main__":
    main()
