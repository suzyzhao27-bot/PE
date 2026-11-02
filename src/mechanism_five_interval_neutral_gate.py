"""Check constant-rate necessity across the five frozen Book2 intervals."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
OUT = ROOT / "data/processed/salt_transport/25_五区间G0恒定速率诊断_20261031.csv"
FIELDS = ("scenario", "sample_id", "lambda_min_s_inv", "lambda_min_window",
          "lambda_max_s_inv", "lambda_max_window", "lambda_max_over_min",
          "delta_min_for_common_rho", "common_rho_at_delta_0p20",
          "scope")


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    if len(lookup) != len(source):
        raise ValueError("Duplicate interval identity")

    output = []
    for scenario in ("C1_T85", "C2_T85", "C1_T80", "C2_T80"):
        for sample in ("B65", "B700", "B2600", "R65", "R700", "R2600"):
            rates = []
            for window in ("I1", "I2", "I3", "I4", "I5"):
                row = lookup[(scenario, window, sample)]
                rate = float(row["lambda_s_inv"])
                if row["status"] != "VALID" or rate <= 0:
                    raise ValueError(f"Invalid interval: {scenario}, {window}, {sample}")
                rates.append((rate, window))
            low, high = min(rates), max(rates)
            delta_min = (high[0] - low[0]) / (high[0] + low[0])
            output.append({
                "scenario": scenario, "sample_id": sample,
                "lambda_min_s_inv": low[0], "lambda_min_window": low[1],
                "lambda_max_s_inv": high[0], "lambda_max_window": high[1],
                "lambda_max_over_min": high[0] / low[0],
                "delta_min_for_common_rho": delta_min,
                "common_rho_at_delta_0p20": delta_min <= .20,
                "scope": "FIVE_INTERVAL_G0_CONSTANT_RATE_NECESSITY_ONLY",
            })
    if len(output) != 24:
        raise ValueError("Expected 24 scenario-curve diagnostics")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"diagnostics={len(output)} no_common_rho_at_delta_0p20="
          f"{sum(r['common_rho_at_delta_0p20'] is False for r in output)}")


if __name__ == "__main__":
    main()
