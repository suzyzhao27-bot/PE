"""Measure conditional two-probe geometry of frozen W/C/WC response contrasts."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
SOURCE = DATA / "22_W_C_WC_处理响应差距_20261028.csv"
OUT = DATA / "22_W_C_WC_双探针方向几何_20261029.csv"
FIELDS = ("h_nm", "salt_feed_mM", "external", "time_min",
          "W_minus_C_salt_uM", "W_minus_C_urea_uM",
          "W_minus_WC_salt_uM", "W_minus_WC_urea_uM",
          "contrast_determinant_uM2", "absolute_angle_sine",
          "core_candidate_if_low_1mM", "core_candidate_if_low_10mM",
          "numeric_status", "experimental_identifiability")


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    if len(source) != 162:
        raise ArithmeticError("Unexpected comparison size")
    groups = defaultdict(dict)
    for row in source:
        key = (int(row["h_nm"]), int(row["salt_feed_mM"]), row["external"], int(row["time_min"]))
        pair = row["mechanism_a"], row["mechanism_b"]
        if pair in groups[key]:
            raise ArithmeticError("Duplicate contrast")
        groups[key][pair] = row
    if len(groups) != 54:
        raise ArithmeticError("Expected 54 condition-time strata")
    output = []
    for (h, feed, external, minute), contrasts in sorted(groups.items()):
        if set(contrasts) != {("W", "C"), ("W", "WC"), ("C", "WC")}:
            raise ArithmeticError("Missing pairwise contrast")
        wc, wwc = contrasts["W", "C"], contrasts["W", "WC"]
        ds1 = float(wc["salt_response_a_uM"]) - float(wc["salt_response_b_uM"])
        du1 = float(wc["urea_response_a_uM"]) - float(wc["urea_response_b_uM"])
        ds2 = float(wwc["salt_response_a_uM"]) - float(wwc["salt_response_b_uM"])
        du2 = float(wwc["urea_response_a_uM"]) - float(wwc["urea_response_b_uM"])
        if abs(du2) > 1e-9:
            raise ArithmeticError("W/WC should differ only in salt for this fixed example")
        determinant = ds1 * du2 - du1 * ds2
        norm_product = math.hypot(ds1, du1) * math.hypot(ds2, du2)
        sine = abs(determinant) / norm_product if norm_product else 0.0
        if not (0 <= sine <= 1 + 1e-12):
            raise ArithmeticError("Invalid contrast angle")
        candidate = h in (65, 2600) or (h == 700 and feed == 600)
        output.append({
            "h_nm": h, "salt_feed_mM": feed, "external": external,
            "time_min": minute, "W_minus_C_salt_uM": ds1,
            "W_minus_C_urea_uM": du1, "W_minus_WC_salt_uM": ds2,
            "W_minus_WC_urea_uM": du2,
            "contrast_determinant_uM2": determinant,
            "absolute_angle_sine": sine,
            "core_candidate_if_low_1mM": int(candidate and feed != 10),
            "core_candidate_if_low_10mM": int(candidate and feed != 1),
            "numeric_status": "NUMERIC_UNRESOLVED_T12" if feed == 10 else "MODEL_ONLY_TARGET_CHECK",
            "experimental_identifiability": "UNDETERMINED_ACTUAL_NOISE_AND_LOQ_NOT_MEASURED",
        })
    if len(output) != 54 or sum(r["core_candidate_if_low_1mM"] for r in output) != 30 or sum(r["core_candidate_if_low_10mM"] for r in output) != 30:
        raise ArithmeticError("Unexpected fixed-core condition-time scope")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    for feed in (1, 10, 600):
        values = [r["absolute_angle_sine"] for r in output if r["salt_feed_mM"] == feed]
        print(f"feed={feed}mM n={len(values)} sine_range={min(values):.6g}..{max(values):.6g}")


if __name__ == "__main__":
    main()
