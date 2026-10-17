"""Plot two audited failure checks for the v5 report; no model fitting."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


THICKNESS = (65, 120, 200, 350, 700, 2600)
COMMON_ROWS = ((0, 0), (1, 20), (2, 40), (3, 60), (5, 100), (6, 120))
CURVE_COLUMNS = {200: "col_5_E", 350: "col_6_F"}


def load_loo(path):
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source))
    if [float(row["held_h_nm"]) for row in rows] != list(THICKNESS):
        raise ValueError("Expected the six frozen held-out thickness rows in order")
    ratios = []
    for row in rows:
        observed = float(row["observed"])
        predicted = float(row["predicted"])
        ratio = float(row["prediction_over_observed"])
        if not all(math.isfinite(value) and value > 0 for value in (observed, predicted, ratio)):
            raise ValueError("Held-out results must be finite and positive")
        if not math.isclose(predicted / observed, ratio, rel_tol=1e-10):
            raise ValueError("Held-out ratio disagrees with predicted/observed")
        ratios.append(ratio)
    return ratios


def load_shapes(path):
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != 7:
        raise ValueError("Expected seven frozen Book1 rows")
    for index, time in COMMON_ROWS:
        if any(float(rows[index][key]) != time for key in ("col_1_A", "col_3_C")):
            raise ValueError(f"A/C time labels disagree at row {index}")
    shapes = {}
    for thickness, column in CURVE_COLUMNS.items():
        start = float(rows[0][column])
        delta = float(rows[6][column]) - start
        if not math.isfinite(delta) or delta <= 0:
            raise ValueError(f"Invalid endpoint increment at {thickness} nm")
        values = [(float(rows[index][column]) - start) / delta for index, _ in COMMON_ROWS]
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"Invalid normalized value at {thickness} nm")
        shapes[thickness] = values
    return shapes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loo-csv", required=True, type=Path)
    parser.add_argument("--book1-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    ratios = load_loo(args.loo_csv)
    shapes = load_shapes(args.book1_csv)
    times = [time for _, time in COMMON_ROWS]

    fig, (left, right) = plt.subplots(1, 2, figsize=(10.8, 4.7))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.24, top=0.84, wspace=0.30)
    left.axhline(1, color="#747b83", linewidth=1.5, linestyle="--", label="Exact held-out match")
    left.scatter(range(6), ratios, color="#17658d", s=52, zorder=3, label="Reference r=0, leave-one-out")
    left.scatter([5], [ratios[5]], color="#c65436", s=72, zorder=4)
    left.annotate(f"{ratios[5]:.1f}×", (5, ratios[5]), xytext=(-8, 9),
                  textcoords="offset points", ha="right", color="#a84228")
    left.set_xticks(range(6), [str(value) for value in THICKNESS])
    left.set_xlabel("Held-out thickness label (nm)")
    left.set_ylabel("Predicted / observed increment")
    left.set_yscale("log")
    left.set_ylim(0.1, 1200)
    left.set_title("A  Held-out endpoint check", loc="left", fontweight="bold")
    left.grid(axis="y", alpha=0.22)
    left.legend(loc="upper left", fontsize=8, frameon=False)

    right.plot(times, [time / 120 for time in times], color="#747b83", linewidth=1.8,
               linestyle="--", label="Constant-P finite-bath lower bound")
    for thickness, color, marker in ((200, "#c65436", "o"), (350, "#17658d", "s")):
        right.plot(times, shapes[thickness], color=color, marker=marker, linewidth=1.8,
                   markersize=5, label=f"Old {thickness} nm curve")
    right.set_xticks(times)
    right.set_xlim(0, 120)
    right.set_ylim(-0.04, 1.08)
    right.set_xlabel("Common time label (min)")
    right.set_ylabel("Normalized conductivity increment")
    right.set_title("B  Conditional time-shape check", loc="left", fontweight="bold")
    right.grid(alpha=0.22)
    right.legend(loc="upper left", fontsize=8, frameon=False)

    fig.text(0.08, 0.93, "Where the shared effective-field model falls short", fontsize=14, fontweight="bold")
    fig.text(0.08, 0.08, "Six old endpoints; 200/350 nm traces use only shared 0, 20, 40, 60, 100, 120 min rows."
             "  No error bars or dynamic-mechanism fit.", fontsize=8, color="#454b52")
    fig.savefig(args.output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
