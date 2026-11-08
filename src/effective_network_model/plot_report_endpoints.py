"""Plot the old thickness series as conditional relative salt-transfer proxies."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


THICKNESS = (65, 120, 200, 350, 700, 2600)


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as source:
        return list(csv.DictReader(source))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint-csv", required=True, type=Path)
    parser.add_argument("--curve-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    endpoints = read_rows(args.endpoint_csv)
    curve = read_rows(args.curve_csv)
    if [float(row["h_nm"]) for row in endpoints] != list(THICKNESS):
        raise ValueError("Expected six frozen thickness endpoints in order")
    observed = [float(row["observed_delta_kappa"]) for row in endpoints]
    predicted = [float(row["predicted_delta_kappa"]) for row in endpoints]
    if not all(math.isfinite(value) and value > 0 for value in observed + predicted):
        raise ValueError("Endpoint increments must be finite and positive")
    for row, obs, pred in zip(endpoints, observed, predicted):
        if not math.isclose(float(row["prediction_over_observed"]), pred / obs, rel_tol=1e-10):
            raise ValueError("Endpoint ratio disagrees with prediction/observation")
    coordinates = [(float(row["h_nm"]), float(row["predicted_delta_kappa"])) for row in curve]
    if len(coordinates) < 2 or any(not math.isfinite(h) or not math.isfinite(y) or h <= 0 or y <= 0
                                   for h, y in coordinates):
        raise ValueError("Fit curve must have at least two finite positive points")
    if any(next_h <= h for (h, _), (next_h, _) in zip(coordinates, coordinates[1:])):
        raise ValueError("Fit curve thickness must strictly increase")
    within = [(h, y) for h, y in coordinates if THICKNESS[0] <= h <= THICKNESS[-1]]
    if len(within) < 2:
        raise ValueError("Fit curve does not span the old thickness-label interval")

    reference_observed = observed[-1]
    reference_predicted = predicted[-1]
    observed_relative = [value / reference_observed for value in observed]
    predicted_relative = [value / reference_predicted for value in predicted]
    curve_relative = [(h, y / reference_predicted) for h, y in within]

    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.83, bottom=0.27)
    ax.plot([h for h, _ in curve_relative], [y for _, y in curve_relative],
            color="#17658d", linewidth=2.2, label="Shared r=0 fit (all six points)")
    ax.scatter(THICKNESS, observed_relative, color="#202a32", s=68, zorder=3,
               label="Old 0–120 min increments / 2600 nm")
    ax.scatter(THICKNESS, predicted_relative, facecolors="white", edgecolors="#17658d",
               s=55, linewidths=1.5, zorder=3, label="Fit at old thickness labels")
    ax.axhline(1, color="#777d82", linewidth=0.8, linestyle=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(60, 2800)
    ax.set_ylim(0.005, 1.6)
    ax.set_xticks(THICKNESS, [str(h) for h in THICKNESS])
    ax.set_xlabel("Old thickness label (nm)")
    ax.set_ylabel("Relative 0–120 min increment (2600 nm = 1)")
    ax.grid(which="major", alpha=0.2)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.text(0.14, 0.92, "Figure 1 | Old thickness trend and relative transfer", fontsize=14,
             fontweight="bold")
    fig.text(0.14, 0.13, "Old points / old 2600 nm; model / its own 2600 nm. Both equal 1 there by definition."
             "\nRelative transfer requires comparable calibration, geometry and volume history; line = six-point training fit.",
             fontsize=8, color="#454b52")
    fig.savefig(args.output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
