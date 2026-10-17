"""Plot only the audited six old endpoints and the shared r=0 training fit."""

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

    fig, ax = plt.subplots(figsize=(8.5, 5.3))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.84, bottom=0.23)
    ax.plot([h for h, _ in within], [y for _, y in within], color="#17658d", linewidth=2.2,
            label="Shared r=0 model, trained on all six endpoints")
    ax.scatter(THICKNESS, observed, color="#202a32", s=65, zorder=3,
               label="Old 0–120 min conductivity increments")
    ax.scatter(THICKNESS, predicted, facecolors="white", edgecolors="#17658d", s=53,
               linewidths=1.5, zorder=3, label="Model at the six labels")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(60, 2800)
    ax.set_ylim(7, 2100)
    ax.set_xticks(THICKNESS, [str(h) for h in THICKNESS])
    ax.set_xlabel("Old thickness label (nm)")
    ax.set_ylabel("0–120 min conductivity increment (inherited table units)")
    ax.grid(which="major", alpha=0.2)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.text(0.13, 0.92, "Old endpoints and the shared effective-field fit", fontsize=14,
             fontweight="bold")
    fig.text(0.13, 0.08, "The line uses all six old endpoints for fitting; only the 65–2600 nm label interval is shown."
             "  No new samples or error bars.", fontsize=8, color="#454b52")
    fig.savefig(args.output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
