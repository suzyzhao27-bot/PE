"""Translate fixed old-device urea envelopes into conditional assay targets."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/processed/salt_transport/25_论文裁剪D尿素浓度条件范围_20261031.csv"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪D检测能力需求_20261031.csv"
FIELDS = ("scenario", "h_nm", "time_min", "B_lower_mM", "B_upper_mM",
          "R_lower_mM", "R_upper_mM", "R_minus_B_gap_mM",
          "predicted_ranges_disjoint", "LOQ_both_below_uM",
          "LOQ_both_for_5x_below_uM", "planned_60min_5x_at_target_1uM",
          "symmetric_absolute_error_below_uM_for_separation", "scope")


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        rows = [r for r in csv.DictReader(stream)
                if r["external_relation"] == "DIFFUSION_SCALED_SHARED"]
    lookup = {(r["scenario"], int(r["h_nm"]), int(r["time_min"]), r["sample_id"][0]): r
              for r in rows}
    if len(rows) != 24 or len(lookup) != 24:
        raise ValueError("D source identity mismatch")
    output = []
    for scenario in ("C1_T85", "C2_T85"):
        for h in (65, 700, 2600):
            for time in (60, 100):
                b = lookup[(scenario, h, time, "B")]
                r = lookup[(scenario, h, time, "R")]
                blo, bhi = float(b["urea_receiver_min_mM"]), float(b["urea_receiver_max_mM"])
                rlo, rhi = float(r["urea_receiver_min_mM"]), float(r["urea_receiver_max_mM"])
                if not (0 < blo <= bhi < 5 and 0 < rlo <= rhi < 5):
                    raise ValueError("Invalid urea prediction interval")
                gap = rlo - bhi
                output.append({
                    "scenario": scenario, "h_nm": h, "time_min": time,
                    "B_lower_mM": blo, "B_upper_mM": bhi,
                    "R_lower_mM": rlo, "R_upper_mM": rhi,
                    "R_minus_B_gap_mM": gap,
                    "predicted_ranges_disjoint": gap > 0,
                    "LOQ_both_below_uM": min(blo, rlo) * 1000,
                    "LOQ_both_for_5x_below_uM": min(blo, rlo) * 200,
                    "planned_60min_5x_at_target_1uM": (
                        min(blo, rlo) * 1000 >= 5 if time == 60 else "OLD_ONLY_NOT_PLANNED_SAMPLE"),
                    "symmetric_absolute_error_below_uM_for_separation": (
                        gap * 500 if gap > 0 else "NO_SEPARATION_GUARANTEE"),
                    "scope": "CONDITIONAL_OLD_DEVICE_ASSAY_TARGET_NOT_VALIDATED_LOQ",
                })
    if len(output) != 12:
        raise ValueError("Expected twelve assay targets")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(output)} disjoint="
          f"{sum(r['predicted_ranges_disjoint'] for r in output)}")


if __name__ == "__main__":
    main()
