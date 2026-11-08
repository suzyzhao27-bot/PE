"""Propagate four independent post-freeze PB points into 10 mM mechanism gaps."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_pb_postfreeze_10mM机制差距影响_20261029.csv"
FIELDS = ("h_nm", "external", "time_min", "mechanism_a", "mechanism_b",
          "frozen_salt_gap_uM", "augmented_salt_gap_uM", "signed_gap_change_uM",
          "absolute_fractional_gap_change", "signed_order_unchanged", "numeric_status")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    curves = read("22_pb_postfreeze_10mM完整时序影响_20261029.csv")
    comparisons = read("22_W_C_WC_处理响应差距_20261028.csv")
    if len(curves) != 48 or len(comparisons) != 162:
        raise ArithmeticError("Unexpected source size")
    change = {}
    for row in curves:
        key = (int(float(row["h_nm"])), row["external"], int(row["time_min"]), row["mechanism"])
        if key in change:
            raise ArithmeticError("Duplicate changed curve state")
        change[key] = float(row["signed_change_uM"])
    output = []
    for row in comparisons:
        if int(row["salt_feed_mM"]) != 10:
            continue
        h, external, minute = int(row["h_nm"]), row["external"], int(row["time_min"])
        a, b = row["mechanism_a"], row["mechanism_b"]
        old_a, old_b = float(row["salt_response_a_uM"]), float(row["salt_response_b_uM"])
        new_a = old_a - (change[h, external, minute, a] if a != "W" else 0.0)
        new_b = old_b - (change[h, external, minute, b] if b != "W" else 0.0)
        old_gap, new_gap = abs(old_a - old_b), abs(new_a - new_b)
        if abs(old_gap - float(row["salt_gap_uM"])) > 1e-8 or old_gap <= 0:
            raise ArithmeticError("Frozen gap did not reproduce")
        output.append({
            "h_nm": h, "external": external, "time_min": minute,
            "mechanism_a": a, "mechanism_b": b,
            "frozen_salt_gap_uM": old_gap, "augmented_salt_gap_uM": new_gap,
            "signed_gap_change_uM": new_gap - old_gap,
            "absolute_fractional_gap_change": abs(new_gap - old_gap) / old_gap,
            "signed_order_unchanged": int((old_a - old_b) * (new_a - new_b) > 0),
            "numeric_status": "POSTFREEZE_SPOTCHECK_T12_UNRESOLVED",
        })
    if len(output) != 54:
        raise ArithmeticError("Expected 54 charged-feed comparisons")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"comparisons=54 order_unchanged={sum(r['signed_order_unchanged'] for r in output)} "
          f"max_absolute_gap_change_uM={max(abs(r['signed_gap_change_uM']) for r in output):.8g} "
          f"max_fractional_gap_change={max(r['absolute_fractional_gap_change'] for r in output):.8g}")


if __name__ == "__main__":
    main()
