"""Propagate 31/16-node trajectory differences into frozen 10 mM salt comparisons."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
CURVES = DATA / "22_T12_10mM完整时序差_20261029.csv"
COMPARISONS = DATA / "22_W_C_WC_处理响应差距_20261028.csv"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = DATA / "22_T12_10mM机制差距敏感性_20261029.csv"
FIELDS = ("h_nm", "external", "time_min", "mechanism_a", "mechanism_b",
          "salt_response_a_31_uM", "salt_response_b_31_uM", "salt_gap_31_uM",
          "salt_response_a_16_uM", "salt_response_b_16_uM", "salt_gap_16_uM",
          "gap_change_uM", "signed_order_unchanged", "numeric_status")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    curves = rows(CURVES)
    comparisons = rows(COMPARISONS)
    if len(curves) != 48 or len(comparisons) != 162:
        raise ArithmeticError("Unexpected input size")
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    charged_pairs = [p for p in pairs if p["c_feed_mM"] == 10 and p["mechanism"] in ("C", "WC")]
    if len(charged_pairs) != 12 or any(p["sigma_B_C_m2"] != -.008 or p["sigma_E_C_m2"] != 0 for p in charged_pairs):
        raise ArithmeticError("B-only charged-state assumption changed")
    adjustment = {}
    for row in curves:
        key = (int(float(row["h_nm"])), row["external"], int(row["time_min"]), row["mechanism"])
        if key in adjustment:
            raise ArithmeticError("Duplicate charged trajectory")
        adjustment[key] = 1000 * (float(row["receiver_31_mM"]) - float(row["receiver_16_mM"]))

    output = []
    for row in comparisons:
        if int(row["salt_feed_mM"]) != 10:
            continue
        h, ext, minute = int(row["h_nm"]), row["external"], int(row["time_min"])
        a, b = row["mechanism_a"], row["mechanism_b"]
        original_a = float(row["salt_response_a_uM"])
        original_b = float(row["salt_response_b_uM"])
        coarse_a = original_a + (adjustment[h, ext, minute, a] if a != "W" else 0.0)
        coarse_b = original_b + (adjustment[h, ext, minute, b] if b != "W" else 0.0)
        fine_gap = abs(original_a - original_b)
        coarse_gap = abs(coarse_a - coarse_b)
        if abs(fine_gap - float(row["salt_gap_uM"])) > 1e-8:
            raise ArithmeticError("Original salt gap did not reproduce")
        output.append({
            "h_nm": h, "external": ext, "time_min": minute,
            "mechanism_a": a, "mechanism_b": b,
            "salt_response_a_31_uM": original_a,
            "salt_response_b_31_uM": original_b,
            "salt_gap_31_uM": fine_gap,
            "salt_response_a_16_uM": coarse_a,
            "salt_response_b_16_uM": coarse_b,
            "salt_gap_16_uM": coarse_gap,
            "gap_change_uM": coarse_gap - fine_gap,
            "signed_order_unchanged": int((original_a - original_b) * (coarse_a - coarse_b) > 0),
            "numeric_status": "NUMERIC_UNRESOLVED_T12",
        })
    if len(output) != 54 or any(r["signed_order_unchanged"] != 1 for r in output):
        raise ArithmeticError("Unexpected 10 mM comparison count or order reversal")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"comparisons=54 max_absolute_gap_change_uM={max(abs(r['gap_change_uM']) for r in output):.8g} "
          f"order_unchanged={sum(r['signed_order_unchanged'] for r in output)}")


if __name__ == "__main__":
    main()
