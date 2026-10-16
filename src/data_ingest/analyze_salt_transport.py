from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = ROOT / "data" / "processed" / "salt_transport"
FIGURE_ROOT = ROOT / "data" / "processed" / "figures"
MPL_ROOT = ROOT / ".mplconfig"
os.environ.setdefault("MPLCONFIGDIR", str(MPL_ROOT))
MPL_ROOT.mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt  # noqa: E402


def load_rows() -> list[dict[str, str]]:
    with (DATA_ROOT / "salt_transport_tidy.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        return list(csv.DictReader(handle))


def summarize(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["curve_id"]].append(row)

    summaries: list[dict[str, object]] = []
    for curve_id, points in sorted(grouped.items()):
        points.sort(key=lambda item: float(item["time"]))
        time = np.asarray([float(point["time"]) for point in points], dtype=float)
        conductivity = np.asarray(
            [float(point["conductivity"]) for point in points], dtype=float
        )
        slope, intercept = np.polyfit(time, conductivity, deg=1)
        predicted = slope * time + intercept
        residual = conductivity - predicted
        ss_res = float(np.sum(residual**2))
        ss_tot = float(np.sum((conductivity - np.mean(conductivity)) ** 2))
        r_squared = float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot
        monotonic_violations = int(np.sum(np.diff(conductivity) < 0))
        first = points[0]
        summaries.append(
            {
                "curve_id": curve_id,
                "source_opju": first["source_opju"],
                "workbook": first["workbook"],
                "worksheet": first["worksheet"],
                "sample_id": first["sample_id"],
                "condition": first["condition"],
                "n_points": len(points),
                "time_min": float(time[0]),
                "time_max": float(time[-1]),
                "conductivity_start": float(conductivity[0]),
                "conductivity_end": float(conductivity[-1]),
                "conductivity_delta": float(conductivity[-1] - conductivity[0]),
                "ols_slope_conductivity_per_min": float(slope),
                "ols_intercept": float(intercept),
                "ols_r_squared_descriptive_only": r_squared,
                "monotonic_violations": monotonic_violations,
                "unit_evidence": first["unit_evidence"],
                "eligible_for_mechanistic_fit": False,
            }
        )
    return summaries


def write_summary(rows: list[dict[str, object]]) -> None:
    path = DATA_ROOT / "salt_transport_curve_summary.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_exact_duplicates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["curve_id"]].append(row)

    signatures: dict[tuple[tuple[float, float], ...], list[str]] = defaultdict(list)
    for curve_id, points in grouped.items():
        signature = tuple(
            sorted(
                (float(point["time"]), float(point["conductivity"]))
                for point in points
            )
        )
        signatures[signature].append(curve_id)

    duplicates: list[dict[str, str]] = []
    for curve_ids in signatures.values():
        if len(curve_ids) < 2:
            continue
        for duplicate in curve_ids[1:]:
            duplicates.append(
                {
                    "canonical_curve_id": curve_ids[0],
                    "duplicate_curve_id": duplicate,
                    "relationship": "exact same time and conductivity values; do not count as an independent replicate",
                }
            )

    path = DATA_ROOT / "exact_duplicate_curves.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["canonical_curve_id", "duplicate_curve_id", "relationship"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(duplicates)
    return duplicates


def plot(rows: list[dict[str, str]]) -> None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["curve_id"]].append(row)

    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    for axis, condition, title in [
        (axes[0], "2h", "Polymer salt project: 2 h tables"),
        (axes[1], "dry_6h", "Polymer salt project: dry/unwetted table"),
    ]:
        for curve_id, points in sorted(grouped.items()):
            if points[0]["condition"] != condition:
                continue
            points.sort(key=lambda item: float(item["time"]))
            time = [float(point["time"]) for point in points]
            conductivity = [float(point["conductivity"]) for point in points]
            label = f'{points[0]["workbook"]}: {points[0]["sample_id"]}'
            axis.plot(time, conductivity, marker="o", markersize=3, linewidth=1.2, label=label)
        axis.set_title(title)
        axis.set_xlabel("Time (min; inherited units flagged in audit)")
        axis.set_ylabel("Conductivity (μS/cm; inherited units flagged in audit)")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=7, ncol=1 if condition == "dry_6h" else 2)

    fig.suptitle(
        "QC-only reconstruction from Origin worksheets — not concentration, flux, or permeability",
        fontsize=11,
    )
    fig.savefig(FIGURE_ROOT / "salt_transport_origin_qc.png", dpi=220)
    plt.close(fig)


def main() -> None:
    rows = load_rows()
    summaries = summarize(rows)
    write_summary(summaries)
    duplicates = write_exact_duplicates(rows)
    plot(rows)
    max_time = max(float(row["time"]) for row in rows if row["condition"] == "dry_6h")
    print(
        {
            "curves": len(summaries),
            "points": len(rows),
            "exact_duplicate_pairs": len(duplicates),
            "dry_6h_max_time_min": max_time,
            "dry_6h_max_time_hours": max_time / 60.0,
        }
    )


if __name__ == "__main__":
    main()
