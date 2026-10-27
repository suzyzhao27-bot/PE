"""Plot saved six-film B holdouts with ND unstable folds visibly marked."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "data/processed/figures"
assert FIG.is_dir()
os.environ.setdefault("MPLCONFIGDIR", str(FIG))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODELS = ("N0", "NS", "ND")
THICKNESSES = (65, 120, 200, 350, 700, 2600)
COLORS = {"N0": "#777777", "NS": "#2767a8", "ND": "#d17925"}


def main() -> None:
    rows = []
    summaries = {}
    for model in MODELS:
        with (ROOT / "WP4" / f"v7_{model}_film_points.csv").open(
                encoding="utf-8-sig", newline="") as stream:
            points = list(csv.DictReader(stream))
        summary = json.loads((ROOT / "WP4" / f"v7_{model}_film_summary.json").read_text(encoding="utf-8"))
        assert len(points) == 42 and len(summary["folds"]) == 6
        summaries[model] = summary
        folds = {int(f["hold_h_nm"]): f for f in summary["folds"]}
        for point in points:
            h = int(point["hold_h_nm"])
            assert point["in_training"] == "False" and h in folds
            rows.append({**point, "fold_status": folds[h]["status"],
                         "model_overall_status": summary["status"]})
    assert len(rows) == 126
    by_key = {}
    for row in rows:
        key = (int(row["hold_h_nm"]), float(row["time_min"]), row["model"])
        assert key not in by_key
        by_key[key] = row
    for h in THICKNESSES:
        reference = [r for r in rows if r["model"] == "N0" and int(r["hold_h_nm"]) == h]
        assert len(reference) == 7
        for base in reference:
            for model in MODELS[1:]:
                other = by_key[h, float(base["time_min"]), model]
                assert (base["curve_id"], base["point_index"],
                        base["observed_increment_uS_cm"]) == (
                        other["curve_id"], other["point_index"],
                        other["observed_increment_uS_cm"])

    source = FIG / "v7_B_film_holdout_panel.csv"
    with source.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 7.1), sharex=True)
    for ax, h in zip(axes.flat, THICKNESSES):
        reference = [r for r in rows if r["model"] == "N0" and int(r["hold_h_nm"]) == h]
        reference.sort(key=lambda r: float(r["time_min"]))
        for model in MODELS:
            points = [by_key[h, float(r["time_min"]), model] for r in reference]
            unstable = points[0]["fold_status"] != "multi_start_stable"
            ax.plot([float(r["time_min"]) for r in points],
                    [float(r["predicted_increment_uS_cm"]) for r in points],
                    color=COLORS[model], linewidth=1.7,
                    linestyle=":" if unstable else ("--" if model == "N0" else "-"),
                    label=model + (" †" if unstable else ""))
        ax.scatter([float(r["time_min"]) for r in reference[1:]],
                   [float(r["observed_increment_uS_cm"]) for r in reference[1:]],
                   s=30, facecolors="white", edgecolors="black",
                   linewidths=1, zorder=5, label="Observed holdout")
        nd_unstable = by_key[h, 0.0, "ND"]["fold_status"] != "multi_start_stable"
        ax.set_title(f"{h} nm" + ("  †" if nd_unstable else ""))
        ax.set_xlim(0, 125)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#e4e4e4", linewidth=0.6)
    for ax in axes[1]:
        ax.set_xlabel("Time (min)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Conductivity increment (µS/cm)")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.suptitle("B series: predict each membrane after fitting the other five", y=0.985)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               ncol=4, frameon=False)
    fig.text(0.5, 0.025,
             "† ND selected fit lacks a second successful start within 1% loss; lines are conditional predictions.",
             ha="center", fontsize=8)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.84, bottom=0.13,
                        wspace=0.19, hspace=0.23)
    for suffix in ("svg", "png"):
        fig.savefig(FIG / f"v7_B_film_holdout_panel.{suffix}", dpi=220,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {source.name} (126 rows) and SVG/PNG; ND unstable folds: 120, 700, 2600 nm")


if __name__ == "__main__":
    main()
