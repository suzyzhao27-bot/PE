"""Show reported zeta trends; pressure ramps are not independent membranes."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "原始资料/赵珈绮/01_研究框架与资料_v2/05_真实数据_v1核查/zeta_summary.csv"
TABLE = ROOT / "data/processed/zeta_trends.csv"
FIGURE = ROOT / "data/processed/figures/zeta_trends.png"


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    if len(source_rows) != 9:
        raise ValueError("Expected nine reported zeta condition summaries")

    rows = []
    for source in source_rows:
        if int(source["ramp_n"]) != 4 or source["independent_membrane_n"] != "unknown":
            raise ValueError("Zeta replicate interpretation changed")
        rows.append({
            "source_id": source["source_id"],
            "draw_ratio_label": int(source["draw_ratio_label"]),
            "pH_mean": float(source["pH_mean"]),
            "zeta_mean_mV": float(source["zeta_mean_mV_reported"]),
            "zeta_sd_four_ramps_ddof0_mV": float(source["zeta_stdev_mV_reported_ddof0"]),
            "ramp_n": 4,
            "independent_membrane_n": "unknown",
            "interpretation": "reported electrokinetic response under original conditions; not charge-source proof",
        })

    fields = list(rows[0])
    with TABLE.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    grouped = defaultdict(list)
    for row in rows:
        grouped[row["draw_ratio_label"]].append(row)
    fig, axis = plt.subplots(figsize=(7.6, 5.2))
    for draw_ratio, points in sorted(grouped.items()):
        points.sort(key=lambda point: point["pH_mean"])
        axis.errorbar(
            [point["pH_mean"] for point in points],
            [point["zeta_mean_mV"] for point in points],
            yerr=[point["zeta_sd_four_ramps_ddof0_mV"] for point in points],
            marker="o", capsize=3, linewidth=1.4, label=f"Draw ratio {draw_ratio}",
        )
    axis.axhline(0, color="0.45", linewidth=0.8)
    axis.set(xlabel="Measured mean pH", ylabel="Reported zeta potential (mV)",
             title="Recorded electrokinetic trends")
    axis.grid(alpha=0.2)
    axis.legend(title="Original labels")
    fig.subplots_adjust(left=0.13, right=0.97, top=0.90, bottom=0.23)
    fig.text(0.5, 0.05,
             "Bars: SD of 4 pressure ramps (ddof=0), not independent membrane SD. Lines guide the eye.",
             ha="center", fontsize=8)
    fig.savefig(FIGURE, dpi=200)
    plt.close(fig)
    print(f"conditions={len(rows)} draw_ratios={len(grouped)} figure={FIGURE.name}")


if __name__ == "__main__":
    main()
