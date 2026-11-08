"""Compare W/C/WC treatment-response gaps across 16/31/35 PB nodes."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_T12_10mM三网格机制差距范围_20261029.csv"
FIELDS = ("h_nm", "external", "time_min", "mechanism_a", "mechanism_b",
          "salt_gap_16_uM", "salt_gap_31_uM", "salt_gap_35_uM",
          "grid_gap_spread_uM", "grid_gap_spread_relative_to_31",
          "order_stable_in_three_grids", "fixed_core_candidate_if_low_10mM",
          "numeric_status")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def key(row: dict[str, str]) -> tuple[int, str, int, str, str]:
    return (int(float(row["h_nm"])), row["external"], int(row["time_min"]),
            row["mechanism_a"], row["mechanism_b"])


def main() -> None:
    coarse = read("22_T12_10mM机制差距敏感性_20261029.csv")
    augmented = read("22_pb_postfreeze_10mM机制差距影响_20261029.csv")
    if len(coarse) != 54 or len(augmented) != 54:
        raise ArithmeticError("Unexpected comparison count")
    refined = {key(r): r for r in augmented}
    if len(refined) != 54:
        raise ArithmeticError("Duplicate augmented comparison")
    output = []
    for row in coarse:
        k = key(row)
        other = refined[k]
        g16, g31, g35 = (float(row["salt_gap_16_uM"]), float(row["salt_gap_31_uM"]),
                         float(other["augmented_salt_gap_uM"]))
        if abs(g31 - float(other["frozen_salt_gap_uM"])) > 1e-8 or g31 <= 0:
            raise ArithmeticError("Frozen gap mismatch")
        spread = max(g16, g31, g35) - min(g16, g31, g35)
        output.append({
            "h_nm": k[0], "external": k[1], "time_min": k[2],
            "mechanism_a": k[3], "mechanism_b": k[4],
            "salt_gap_16_uM": g16, "salt_gap_31_uM": g31,
            "salt_gap_35_uM": g35, "grid_gap_spread_uM": spread,
            "grid_gap_spread_relative_to_31": spread / g31,
            "order_stable_in_three_grids": int(row["signed_order_unchanged"] == "1"
                                              and other["signed_order_unchanged"] == "1"),
            "fixed_core_candidate_if_low_10mM": int(k[0] in (65, 2600)),
            "numeric_status": "INTERPOLATION_GRID_SPREAD_T12_UNRESOLVED",
        })
    if len(output) != 54 or sum(r["fixed_core_candidate_if_low_10mM"] for r in output) != 36:
        raise ArithmeticError("Unexpected fixed-core comparison scope")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"comparisons=54 order_stable={sum(r['order_stable_in_three_grids'] for r in output)} "
          f"max_spread_uM={max(r['grid_gap_spread_uM'] for r in output):.8g} "
          f"max_relative={max(r['grid_gap_spread_relative_to_31'] for r in output):.8g}")


if __name__ == "__main__":
    main()
