"""Compare frozen W1/W2 interval rates without fitting a transport model."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
OUT = ROOT / "data/processed/salt_transport/25_W1W2速率区间诊断_20261030.csv"
SAMPLES = ("B65", "B700", "B2600", "R65", "R700", "R2600")
FIELDS = ("scenario", "sample_id", "lambda_W1_s_inv", "lambda_W2_s_inv",
          "lambda_W2_over_W1", "delta_min_for_rho_overlap",
          "rho_intervals_overlap_at_delta_0p20", "scope")


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in rows}
    if len(lookup) != len(rows):
        raise ValueError("Duplicate frozen rate identity")

    output = []
    for conversion in ("C1", "C2"):
        for sample in SAMPLES:
            rates = []
            for window in ("W1", "W2"):
                keys = [(f"{conversion}_T{time}", window, sample) for time in (85, 80)]
                a, b = (lookup[key] for key in keys)
                if a["status"] != "VALID" or b["status"] != "VALID":
                    raise ValueError(f"Invalid frozen rate: {keys}")
                if a["lambda_s_inv"] != b["lambda_s_inv"]:
                    raise ValueError(f"T80/T85 window rates differ: {keys}")
                rate = float(a["lambda_s_inv"])
                if rate <= 0:
                    raise ValueError(f"Nonpositive rate: {keys[0]}")
                rates.append(rate)
            first, second = rates
            threshold = abs(second - first) / (second + first)
            output.append({
                "scenario": conversion, "sample_id": sample,
                "lambda_W1_s_inv": first, "lambda_W2_s_inv": second,
                "lambda_W2_over_W1": second / first,
                "delta_min_for_rho_overlap": threshold,
                "rho_intervals_overlap_at_delta_0p20": threshold <= 0.20,
                "scope": "SINGLE_CURVE_INTERVAL_RATE_ONLY",
            })
    if len(output) != 12:
        raise ValueError("Expected twelve rate diagnostics")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"diagnostics={len(output)} nonoverlap_at_delta_0p20="
          f"{sum(r['rho_intervals_overlap_at_delta_0p20'] is False for r in output)}")


if __name__ == "__main__":
    main()
